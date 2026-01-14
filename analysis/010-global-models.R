# EME Linux: 2.2 TB RAM, 2 Intel Xeon Platinum 8462Y+ processors, 64 cores
library('dplyr')     # for data wrangling
library('purrr')     # for functional programing
library('furrr')     # for parallelized functional programing
library('lubridate') # for working with dates
library('mgcv')      # for Genralized Additive Models
library('gratia')    # for plotting GAMs
library('qs')        # for saving R objects quicker than with saveRDS
library('sf')        # for splitting predictions across groups
library('terra')     # for working with rasters
source('analysis/figures/000-default-ggplot-theme.R')

NCORES <- future::availableCores(logical = FALSE) - 2
t_res <- 4
s_res <- 4

#' *only ran neotropic-nearctic (4) for now*
GROUP <-
  st_read('data/ecoregions/groups-polygons.shp', quiet = TRUE) %>%
  filter(group != 'Antarctic') %>%
  pull(group) %>%
  unique() %>%
  tolower() %>%
  gsub(', ', '-', .) %>%
  nth(3) #' *change this as one of 1, 2, 3, or 4*
GROUP
DATE <- paste(Sys.Date()) # make a character to avoid bad conversions
DATE

# fit the mean models for each group ----
# fitting times
#' neotropic-nearctic: 2 days; 1.7 h w/o `ti(year, doy, x, y)` 
#' africa: 
#' palearctic: 
#' indo-malay-oceania-australasia: 
# import data ----
d <- paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
            GROUP, '-t-', t_res, '-s-', s_res, '-ndvi-data.qs') %>%
  qread(nthreads = NCORES)

#' `k` values for spatial smooths
k_s <-
  case_when(GROUP == 'africa' ~ 2000,
            GROUP == 'indo-malay-oceania-australasia' ~ 900,
            GROUP == 'neotropic-nearctic' ~ 1300,
            GROUP == 'palearctic' ~ 2000)

# fit mean model ----
if(length(list.files('H:/GitHub/ndvi-stochasticity/models/global-models',
                     paste0('bam-mean-ndvi-', GROUP))) > 0) {
  m_mu <- qread(list.files(path = 'models/global-models',
                           pattern = paste0('bam-mean-ndvi-', GROUP),
                           full.names = TRUE),
                nthreads = NCORES)
} else {
  gc() # clean up before starting
  
  m_mu <- bam(
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
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    nthreads = NCORES,
    control = gam.control(trace = TRUE),
    samfrac = 0.01) # initial guess using only 1% of data (default = 100%)
  
  qsave(m_mu, paste0('models/global-models/bam-mean-ndvi-', GROUP, '-',
                     DATE, '.qs'), nthreads = NCORES)
}

# find pixel-level variance ----
#' preds for neotropic-nearctic ran with `discrete = FALSE` (34-41 days)
#' preds for all other groups ran with `discrete = TRUE` due to time
#' constraints for the PhD
d$mu_hat <- predict.bam(m_mu, discrete = GROUP != 'neotropic-nearctic',
                        block.size = 1e6, newdata.guaranteed = TRUE) %>%
  as.numeric()

d <- mutate(d, e = (ndvi_aggr - mu_hat), e2 = e^2)

qsave(d,
      paste0('models/global-models/neotropic-nearcric-preds-w-resid-', Sys.Date(),'.qs'),
      nthreads = NCORES)

#' check pixel-level residuals: `dicrete = TRUE` is too coarse but faster
if(FALSE) {
  z <- d %>%
    slice(1:2e5) %>%
    filter(central_date == central_date[1e5])
  
  z %>%
    ggplot(aes(x, y, fill = mu_hat)) +
    geom_raster() +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  z %>%
    mutate(mu_hat = predict.bam(m_mu, discrete = TRUE, newdata = z,
                                type = 'response')) %>%
    ggplot(aes(x, y, fill = mu_hat)) +
    geom_raster() +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  rm(z)
}

# model variance ----
gc() # clean up before starting

# neotropic-nearctic ran in 38 hours
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
  data = d,
  method = 'fREML',
  knots = list(doy = c(0.5, 366.5)),
  discrete = TRUE,
  nthreads = NCORES,
  control = gam.control(trace = TRUE),
  samfrac = 0.01) # initial guess using only 1% of data (default = 100%)

paste0('models/global-models/bam-variance-ndvi-', GROUP,'-',DATE,'.qs') %>%
  qsave(m_s2, ., nthreads = NCORES)

#' *DID NOT RUN: may have bugs*
if(FALSE) {
  # save model summary ----
  sink(file = 'models/global-models/global-model-summaries.txt',
       append = TRUE)
  
  cat('region: ', GROUP, '; parameter: mean\n', sep = '')
  print(summary(m_mu))
  
  cat('\n- - - - -\n\n', sep = '')
  
  cat('region: ', GROUP, '; parameter: variance\n', sep = '')
  print(summary(m_s2))
  
  cat('\n= = = = = = = = = =\n', sep = '')
  cat('- - - - - - - - - -\n', sep = '')
  cat('= = = = = = = = = =\n\n', sep = '')
  
  sink() # close connection to text file
  
  # save a plot of the smooth terms for the mean model ----
  png(paste0('figures/global-models/mean-ndvi-terms-', GROUP,
             '-t-4-s-4-', DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(m_mu, pages = 1, too.far = 0.05, scale = 0,
       scheme = c(1, 5, 1, 1, 1, 3, 5, 5))
  dev.off()
  
  # save a plot of the smooth terms for the variance model ----
  png(paste0('figures/global-models/variance-ndvi-terms-', g, '-',
             DATE, '.png'),
      width = 12, height = 9, units = 'in', bg = 'white', res = 300)
  plot(m_s2, pages = 1, too.far = 0.05, scale = 0,
       scheme = c(1, 5, 1, 1, 1, 3, 5, 5))
  dev.off()
  
  # write rasters of cleaned NDVI, mean NDVI, and residuals ----
  paste0('data/output/', GROUP, '/',
         c('aggregated-unmodeled-ndvi',
           'estimated-mean-ndvi-fitted-only',
           'squared-ndvi-residuals')) %>%
    map_chr(function(.dir) {
      dir.create(.dir)
      return(paste(.dir, 'created.'))
    })
  
  plan('multisession', workers = NCORES)
  
  d %>%
    select(central_date, x, y, ndvi_aggr, mu_hat, e2) %>%
    nest(r = ! central_date) %>%
    future_map2(as.character(central_date), r, function(.date, .r) {
      .r %>%
        rast() %>%
        terra::`crs<-`('EPSG:4326') %>%
        writeRaster(
          paste0('data/outputs/neotropic-nearctic/',
                 c('aggregated-unmodeled-ndvi/observed-ndvi-',
                   'estimated-mean-ndvi-fitted-only/fitted-',
                   'squared-ndvi-residuals/e2-',
                   'DENVar/denvar-')),
          .date, '.tif'))
  return(paste('Exported', .date))
    }, .options = furrr_options(seed = NULL), .progress = TRUE)

}
