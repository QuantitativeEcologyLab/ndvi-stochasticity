# EME Linux: 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
library('dplyr')   # for data wrangling
library('mgcv')    # for Genralized Additive Models
library('gratia')  # for plotting GAMs
library('qs')      # for saving R objects quicker than with saveRDS
source('analysis/figures/000-default-ggplot-theme.R')

NCORES <- future::availableCores(logical = FALSE) - 2

GROUPS <-
  sf::st_read('data/ecoregions/ecoregions-polygons.shp', quiet = TRUE) %>%
  filter(group != 'Antarctic') %>%
  pull(group) %>%
  unique() %>%
  tolower() %>%
  gsub(', ', '-', .)
GROUPS

# set up the data for the model
map(GROUPS, function(GROUP) {
  gc() # clean up before starting
  
  d <- qread(paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
                    GROUP, '-t-2-s-2-ndvi-data.qs'), nthreads = NCORES)
  
  DATE <- paste(Sys.Date()) # make a character to avoid bad conversions
  
  k_s <-
    case_when(GROUP == 'africa' ~ 2000,
              GROUP == 'indo-malay-oceania-australasia' ~ 900,
              GROUP == 'neotropic-nearctic' ~ 1300,
              GROUP == 'palearctic' ~ 2000)
  
  m_mu <- bam(
    ndvi ~
      # marginal smooths
      s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
      s(y, x, bs = 'sos', k = k_s) + # smooth of space
      s(year, bs = 'cr', k = 20) + # year effect
      s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
      s(elev_m, bs = 'cr', k = 10) + # elevation effect
      # tensor product interaction smooths
      ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
      ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(10, k_s / 2)) +
      ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 2)),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    nthreads = NCORES,
    control = gam.control(trace = TRUE),
    samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
  
  qsave(m_mu, paste0('models/global-models/bam-mean-ndvi-', GROUP, '-',
                     DATE, '.qs'), nthreads = NCORES)
  
  sink(file = 'models/global-models/global-model-summaries.txt',
       append = TRUE)
  cat('region: ', GROUP, '; parameter: mean; date: ', DATE, '\n',
      sep = '')
  summary(m_mu)
  cat('\n= = = = = = = = = = = = = = = = = = = = = = = = =\n\n', sep = '')
  sink() # close connection to text file
  
  png(paste0('figures/global-models/mean-ndvi-terms-', GROUP, '-',
             DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(m_mu, pages = 1, too.far = 0.01, scheme = c(1, 5, 1, 1, 1, 3, 5, 5),
       scale = 0)
  dev.off()
  
  # find pixel-level variance ----
  # model has an identity link function: link = response
  d <- bind_cols(d,
                 predict.bam(m_mu, newdata = d, type = 'response',
                             se.fit = TRUE, newdata.guaranteed = TRUE,
                             nthreads = NCORES, block.size = 5e6) %>%
                   as.data.frame() %>%
                   rename(mu_hat = fit,
                          se_mu_hat = se.fit)) %>%
    mutate(e_2 = (ndvi_aggr - mu_hat)^2)
  
  qsave(d, paste0('data/bam-var-ndvi-data-', GROUP, '-', DATE, '.qs'),
        nthreads = NCORES)
  
  # check pixel-level mean residuals
  if(FALSE) {
    d %>%
      rename(long = x, lat = y) %>%
      group_by(long, lat) %>%
      summarize(mean_e = mean(e, na.rm = TRUE), .groups = 'drop') %>%
      pull(mean_e) %>%
      hist(breaks = 100)
    
    d %>%
      rename(long = x, lat = y) %>%
      group_by(long, lat) %>%
      summarize(median_e = median(e, na.rm = TRUE), .groups = 'drop') %>%
      pull(median_e) %>%
      hist(breaks = 100)
  }
  
  # model variance
  m_s2 <- bam(
    e2 ~
      # marginal smooths
      s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
      s(y, x, bs = 'sos', k = k_s) + # smooth of space
      s(year, bs = 'cr', k = 20) + # year effect
      s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
      s(elev_m, bs = 'cr', k = 10) + # elevation effect
      # tensor product interaction smooths
      ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
      ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(10, k_s / 2)) +
      ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 2)),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    nthreads = NCORES,
    control = gam.control(trace = TRUE),
    samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
  
  qsave(m_s2, paste0('models/global-models/bam-variance-ndvi-', GROUP,
                     '-', DATE, '.qs'), nthreads = NCORES)
  
  sink(file = 'models/global-models/global-model-summaries.txt',
       append = TRUE)
  cat('region: ', GROUP, '; parameter: variance; date: ', DATE, '\n',
      sep = '')
  summary(m_s2)
  cat('\n= = = = = = = = = = = = = = = = = = = = = = = = =\n\n', sep = '')
  sink() # close connection to text file
  
  png(paste0('figures/global-models/variance-ndvi-terms-', GROUP, '-',
             DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(m_s2, pages = 1, too.far = 0.01, scheme = c(1, 5, 1, 1, 1, 3, 5, 5),
       scale = 0)
  dev.off()
})
