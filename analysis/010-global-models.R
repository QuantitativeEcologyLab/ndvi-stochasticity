library('dplyr')   # for data wrangling
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('gratia')  # for plotting GAMs
source('analysis/figures/000-default-ggplot-theme.R')

eco <- st_read('data/ecoregions/ecoregions-polygons.shp', quiet = TRUE)

GROUPS <-
  unique(eco$group) %>%
  tolower() %>%
  gsub(', ', '-', .)

# set up the data for the model
map(GROUPS, function(GROUP) {
  gc() # clean up before starting
  
  d <- readRDS(paste0('data/ndvi-data-', GROUP, '.rds'))
  
  DATE <- paste(Sys.Date()) # make a character to avoid bad conversions
  
  k_s <-
    case_when(GROUP == 'australasia-indo-malay-oceania' ~ 1e3,
              GROUP == 'antarctic' ~ 1e3,
              GROUP == 'afrotropic' ~ 1e3,
              GROUP == 'nearctic-neotropic' ~ 1e3,
              GROUP == 'palearctic' ~ 1e3)
  
  m_mu <- bam(
    ndvi ~
      # marginal smooths
      s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
      s(y, x, bs = 'sos', k = k_s) + # smooth of space
      s(year, bs = 'cr', k = 20) + # year effect
      s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
      s(elev_m, bs = 'cr', k = 10) + # elevation effect
      # tensor interaction smooths
      ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
      ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(6, k_s / 2)) +
      ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 2)),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    nthreads = 60,
    control = gam.control(trace = TRUE),
    samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
  
  saveRDS(m_mu, paste0('models/global-models/bam-mean-ndvi-', GROUP, '-',
                       DATE, '.rds'))
  
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
  
  # plot each smooth separately
  p_bam_doy <- draw(m_mu, rug = FALSE,
                    select = which(grepl('doy', smooths(m_mu))))
  ggsave(paste0('figures/bam-mean-ndvi-doy-', GROUP, '-', DATE, '.png'),
         plot = p_bam_doy, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
  
  p_bam_y <- draw(m_mu, rug = FALSE,
                  select = which(grepl('year', smooths(m_mu))))
  ggsave(paste0('figures/bam-mean-ndvi-year-', GROUP, '-', DATE, '.png'),
         plot = p_bam_y, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
  
  p_bam_elev <- draw(m_mu, rug = FALSE,
                     select = which(smooths(m_mu) == 's(elevation_m)'))
  ggsave(paste0('figures/bam-mean-ndvi-elev_m-', GROUP, '-', DATE, '.png'),
         plot = p_bam_elev, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
  
  # find pixel-level variance ----
  # model has an link function: link = response
  d$mu_hat <- fitted(m_mu)
  d$e_2 <- residuals(m_mu)^2
  d$e_2 <- d$e^2 # Var(Y) = E((Y - E(Y))^2) = E(e^2)
  saveRDS(d, paste0('data/bam-var-ndvi-data-', GROUP, '-', DATE, '.rds'))
  
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
      # tensor interaction smooths
      ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
      ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(6, k_s / 2)) +
      ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 2)),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    nthreads = 60,
    control = gam.control(trace = TRUE),
    samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
  
  saveRDS(m_s2, paste0('models/global-models/bam-variance-ndvi-', GROUP,
                       '-', DATE, '.rds'))
  
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
  
  # plot each smooth separately
  p_bam_doy <- draw(m_s2, rug = FALSE,
                    select = which(grepl('doy', smooths(m_s2))))
  ggsave(paste0('figures/bam-variance-ndvi-doy-', GROUP, '-', DATE, '.png'),
         plot = p_bam_doy, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
  
  p_bam_y <- draw(m_s2, rug = FALSE,
                  select = which(grepl('year', smooths(m_s2))))
  ggsave(paste0('figures/bam-variance-ndvi-year-', GROUP, '-', DATE, '.png'),
         plot = p_bam_y, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
  
  p_bam_elev <- draw(m_s2, rug = FALSE,
                     select = which(smooths(m_s2) == 's(elevation_m)'))
  ggsave(paste0('figures/bam-variance-ndvi-elev_m-', GROUP, '-', DATE, '.png'),
         plot = p_bam_elev, width = 8, height = 6, units = 'in', dpi = 300,
         bg = 'white')
})
