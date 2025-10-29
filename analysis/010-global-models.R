# EME Linux: 2.2 TB RAM, 2 Intel Xeon Platinum 8462Y+ processors, 64 cores
library('dplyr')     # for data wrangling
library('purrr')     # for functional programing
library('lubridate') # for working with dates
library('mgcv')      # for Genralized Additive Models
library('gratia')    # for plotting GAMs
library('qs')        # for saving R objects quicker than with saveRDS
source('analysis/figures/000-default-ggplot-theme.R')

NCORES <- future::availableCores(logical = FALSE) - 2
t_res <- 4
s_res <- 4

GROUPS <-
  sf::st_read('data/ecoregions/groups-polygons.shp', quiet = TRUE) %>%
  filter(group != 'Antarctic') %>%
  pull(group) %>%
  unique() %>%
  tolower() %>%
  gsub(', ', '-', .)
GROUPS <- GROUPS[4:1]
GROUPS
DATE <- paste(Sys.Date()) # make a character to avoid bad conversions
DATE

# fit the mean models for each group ----
# 3x aggregation maxes out RAM for neotropic-nearctic model
# if 4x crashes, reduce k of tis to k_s/4
# fitting times
# neotropic-nearctic: 21-min setup, 22-min gam.setup, 1 h fit; ~1.7 h total
# africa: 
# palearctic: 
# islands: 
models <-
  tibble(group = GROUPS,
         data = map(group, function(.g) {
           paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
                  .g, '-t-', t_res, '-s-', s_res, '-ndvi-data.qs') %>%
             qread(nthreads = NCORES)
         }),
         m_mu = map2(group, data, function(.g, .d) {
           gc() # clean up before starting
           
           k_s <-
             case_when(.g == 'africa' ~ 2000,
                       .g == 'indo-malay-oceania-australasia' ~ 900,
                       .g == 'neotropic-nearctic' ~ 1300,
                       .g == 'palearctic' ~ 2000)
           
           .m_mu <- bam(
             ndvi_aggr ~
               # marginal smooths
               s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
               s(y, x, bs = 'sos', k = k_s) + # smooth of space
               s(year, bs = 'cr', k = 20) + # year effect
               s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
               s(elev_m, bs = 'cr', k = 10) + # elevation effect
               # tensor product interaction smooths
               # seasonal trends vary over the years
               ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
               # yearly trends vary spatially
               ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(10, k_s / 4)) +
               # seasonal trends vary spatially
               ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 4)) +
               # change in seasonal trends over the years varies spatially
               # (e.g., seasonal trends in the arctic change differently
               # from those at the equator)
               ti(year, doy, y, x, bs = c('cr', 'cc', 'sos'), d = c(1, 1, 2),
                  k = c(10, 10, k_s / 4)),
             family = gaussian(),
             data = .d,
             method = 'fREML',
             knots = list(doy = c(0.5, 366.5)),
             discrete = TRUE,
             nthreads = NCORES,
             control = gam.control(trace = TRUE),
             samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
           
           qsave(.m_mu, paste0('models/global-models/bam-mean-ndvi-', .g, '-',
                               DATE, '.qs'), nthreads = NCORES)
           return(.m_mu)
         }, .progress = TRUE))

# find pixel-level variance ----
# model has an identity link function: link = response
dates <- readRDS('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds') %>%
  mutate(julian = julian(date),
         date_group = as.integer(julian - (julian %% t_res))) %>%
  # dates and number of rasters for each date_group
  group_by(date_group) %>%
  mutate(n_rasters = sum(! is.na(file_name)),
         start_date = min(date),
         central_date = mean(date),
         end_date = max(date),
         doy = yday(central_date) %>% as.integer(),
         year = year(central_date) %>% as.integer()) %>%
  ungroup()



# d <- mutate(d,
#             mu_hat = predict(m_mu, newdata = d, type = 'response',
#                              discrete = FALSE, se.fit = FALSE,
#                              block.size = 1e6),
#             e2 = (ndvi_aggr - mu_hat)^2)
# 
# qsave(d, paste0('data/bam-var-ndvi-data-', .g, '-t-4-s-4-', DATE, '.qs'),
#       nthreads = NCORES)

# check pixel-level residuals
if(FALSE) {
  z <- d %>%
    slice(1:2e5) %>%
    filter(central_date == central_date[1e5])
  
  z %>%
    ggplot(aes(x, y, fill = mu_hat)) +
    geom_raster() +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  z %>%
    mutate(mu_hat = predict.bam(m_mu, discrete = FALSE, newdata = z,
                                type = 'response')) %>%
    ggplot(aes(x, y, fill = mu_hat)) +
    geom_raster() +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  rm(z)
}

# model variance
models <- models %>%
  mutate(m_s2 = map2(group, data, function(.g, .d) {
    gc() # clean up before starting
    
    m_s2 <- bam(
      e2 ~
        # marginal smooths
        s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
        s(y, x, bs = 'sos', k = k_s) + # smooth of space
        s(year, bs = 'cr', k = 20) + # year effect
        s(doy, bs = 'cc', k = 10) + # seasonal/day of year effect
        s(elev_m, bs = 'cr', k = 10) + # elevation effect
        s(mu_hat, bs = 'cr', k = 5) + # mean and variance are correlated
        # tensor product interaction smooths
        # tensor product interaction smooths
        # seasonal trends vary over the years
        ti(year, doy, bs = c('cr', 'cc'), k = c(10, 10)) +
        # yearly trends vary spatially
        ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(10, k_s / 4)) +
        # seasonal trends vary spatially
        ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, k_s / 4)) +
        # change in seasonal trends over the years varies spatially
        # (e.g., seasonal trends in the arctic change differently
        # from those at the equator)
        ti(year, doy, y, x, bs = c('cr', 'cc', 'sos'), d = c(1, 1, 2),
           k = c(10, 10, k_s / 4)),
      family = gaussian(),
      data = .d,
      method = 'fREML',
      knots = list(doy = c(0.5, 366.5)),
      discrete = TRUE,
      nthreads = NCORES,
      control = gam.control(trace = TRUE),
      samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
    
    qsave(m_s2, paste0('models/global-models/bam-variance-ndvi-', g,
                       '-', DATE, '.qs'), nthreads = NCORES)
    return(m_s2)
  }))

# save the model summaries and plots
imap(models$group, function(.i, .g) {
  # model summaries
  sink(file = 'models/global-models/global-model-summaries.txt',
       append = TRUE)
  
  cat('region: ', .g, '; parameter: mean\n', sep = '')
  print(summary(models$m_mu[[i]]))
  
  cat('\n- - - - -\n\n', sep = '')
  
  cat('region: ', g, '; parameter: variance\n', sep = '')
  print(summary(models$m_s2[[i]]))
  
  cat('\n= = = = = = = = = =\n', sep = '')
  cat('- - - - - - - - - -\n', sep = '')
  cat('= = = = = = = = = =\n\n', sep = '')
  
  sink() # close connection to text file
  
  # save a plot of the smooth terms
  # mean
  png(paste0('figures/global-models/mean-ndvi-terms-', .g,
             '-t-4-s-4-', DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(models$m_mu[[.i]], pages = 1, too.far = 0.05, scale = 0,
       scheme = c(1, 5, 1, 1, 1, 3, 5, 5))
  dev.off()
  
  # variance
  png(paste0('figures/global-models/variance-ndvi-terms-', g, '-',
             DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(models$m_s2[[.i]], pages = 1, too.far = 0.05, scale = 0,
       scheme = c(1, 5, 1, 1, 1, 3, 5, 5))
  dev.off()
})
