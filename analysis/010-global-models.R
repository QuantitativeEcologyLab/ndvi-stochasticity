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

GROUPS <-
  st_read('data/ecoregions/groups-polygons.shp', quiet = TRUE) %>%
  filter(group != 'Antarctic') %>%
  pull(group) %>%
  unique() %>%
  tolower() %>%
  gsub(', ', '-', .)
GROUPS <- GROUPS[4]#:1]
GROUPS
DATE <- paste(Sys.Date()) # make a character to avoid bad conversions
DATE

# fit the mean models for each group ----
# 3x aggregation maxes out RAM for neotropic-nearctic model
# if 4x crashes, reduce k of tis to k_s/4
# fitting times
#' neotropic-nearctic: 2 days; 1.7 h w/o `ti(year, doy, x, y)` 
#' africa: 
#' palearctic: 
#' islands: 
models <-
  tibble(group = GROUPS,
         data = map(group, function(.g) {
           paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
                  .g, '-t-', t_res, '-s-', s_res, '-ndvi-data.qs') %>%
             qread(nthreads = NCORES)
         }),
         m_mu = map2(group, data, function(.g, .d) {
           if(length(list.files('H:/GitHub/ndvi-stochasticity/models/global-models',
                                paste0('bam-mean-ndvi-', .g))) > 0) {
             .m_mu <-
               qread(list.files('models/global-models',
                                paste0('bam-mean-ndvi-', .g),
                                full.names = TRUE),
                     nthreads = NCORES)
           } else {
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
           }
           return(.m_mu)
         }, .progress = TRUE))

# find pixel-level variance ----
m_mu_nn <- models$m_mu[[1]]
d_nn <- models$data[[1]]
d_nn$mu_hat <- predict.bam(m_mu_nn, discrete = FALSE, block.size = 1e6,
                           newdata.guaranteed = TRUE)
d_nn <- mutate(d_nn, e2 = (ndvi_aggr - mu_hat)^2)

paste0('data/output/neotropic-nearctic/',
       c('aggregated-unmodeled-ndvi',
         'estimated-mean-ndvi-fitted-only',
         'squared-ndvi-residuals')) %>%
  map_chr(function(.dir) {
    dir.create(.dir)
    return(paste(.dir, 'created.'))
  })

plan('multisession', workers = NCORES)

d_nn %>%
  select(central_date, x, y, ndvi_aggr, mu_hat, e2) %>%
  nest(r = ! central_date) %>%
  future_map2(as.character(central_date), r, function(.date, .r) {
    .r %>%
      rast() %>%
      terra::`crs<-`('EPSG:4326')%>%
      writeRaster(
        paste0('data/outputs/neotropic-nearctic/',
               c('aggregated-unmodeled-ndvi/observed-ndvi-',
                 'estimated-mean-ndvi-fitted-only/fitted-',
                 'squared-ndvi-residuals/e2-'),
               .date, '.tif'))
    return(paste('Exported', .date))
  }, .options = furrr_options(seed = NULL), .progress = TRUE)

#' *DELETE? CHECK BEFORE RUNNING*
# future_map(raw_file_names, function(.fn) {
#   preds <- rast(.fn) %>%
#     as.data.frame(xy = TRUE, na.rm = TRUE) %>%
#     # add model covariates
#     mutate(.,
#            date =
#              gsub('*.unmodeled-ndvi-', '', .fn) %>%
#              gsub('-t-4-s-4.tif', '', .) %>%
#              as_date(),
#            year = year(date),
#            doy = yday(date),
#            elev_m = extract(r_elev, select(., x, y))[, 2],
#            prop_water = 0,
#            group = select(., x, y) %>%
#              st_as_sf(coords = c('x', 'y')) %>%
#              st_set_crs('EPSG:4326') %>%
#              st_transform(st_crs(group_shp)) %>%
#              mutate(.,
#                     group = st_intersects(., group_shp),
#                     group = as.numeric(group),
#                     group = group_shp$group[group]) %>%
#              pull(group)) %>%
#     nest(group_data = ! group) %>%
#     # add predictions and squared residuals
#     mutate(group_data = map2(group, group_data, function(..g, ..data) {
#       ..model <- models$model[[which(models$group == ..g)]]
#       
#       
#       mutate(..data,
#              mu_hat = predict(m_mu, newdata = d, type = 'response',
#                               discrete = FALSE),
#              e = ndvi_aggr - mu_hat,
#              e2 = (e)^2)
#     })) %>%
#     unnest(group_data)
#   
#   preds %>%
#     select(x, y, mu_hat, e, e2) %>%
#     rast() %>%
#     writeRaster(filename = c(
#       gsub('cleaned-and-aggregated/unmodeled-ndvi',
#            'estimated-means/mu-hat', .fn),
#       gsub('cleaned-and-aggregated/unmodeled-ndvi',
#            'residuals/residuals', .fn),
#       gsub('cleaned-and-aggregated/unmodeled-ndvi',
#            'squared-residuals/squared-residuals', .fn)))
# })

#' # check pixel-level residuals *EDIT OR REMOVE*
# if(FALSE) {
#   z <- d %>%
#     slice(1:2e5) %>%
#     filter(central_date == central_date[1e5])
#   
#   z %>%
#     ggplot(aes(x, y, fill = mu_hat)) +
#     geom_raster() +
#     scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
#   
#   z %>%
#     mutate(mu_hat = predict.bam(m_mu, discrete = FALSE, newdata = z,
#                                 type = 'response')) %>%
#     ggplot(aes(x, y, fill = mu_hat)) +
#     geom_raster() +
#     scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
#   
#   rm(z)
# }

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
    
    qsave(m_s2,
          paste0('models/global-models/bam-variance-ndvi-', g, '-', DATE,'.qs'),
          nthreads = NCORES)
    return(m_s2)
  }))

s2_file_names <- gsub('unmodeled-ndvi', 'DENVar', raw_file_names)
#' *ADD CODE TO SAVE DENVar ESTIMATES*

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
