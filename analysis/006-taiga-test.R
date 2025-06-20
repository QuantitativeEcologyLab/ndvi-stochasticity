# specs of the EME Linux:
# 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
#' some code wrappted in `if(FALSE){...}` because it takes too long to run
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('elevatr')   # for digital elevation models
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('mgcv')      # for GAMs
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
library('gratia')    # for fancy plots of GAMs
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('functions/decode_qa.R')    # for cleaning rasters
source('functions/bit_to_int.R')   # for cleaning rasters
source('functions/is_flagged.R') # for cleaning rasters
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R') # get_legend() from cowplot v. 1.1.3 fails
source('functions/nbs_from_rast.R') # gives a list of neighboring cells

# pick a northern polygon
na <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(WWF_REALM == 'NA') %>%
  st_transform(crs(rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc')))

if(FALSE) {
  # searching for a specific polygon
  # second digit is a 2
  ggplot(st_geometry(na), aes(fill = factor(1:nrow(na) %% 10))) +
    geom_sf() +
    scale_fill_discreterainbow()
  
  # first digit is 7
  ggplot(st_geometry(na), aes(fill = factor(floor(1:nrow(na) / 10)))) +
    geom_sf() +
    scale_fill_discreterainbow()
}

taiga <- na %>%
  slice(72) %>%
  st_cast('POLYGON', warn = FALSE) %>%
  st_geometry() %>%
  st_as_sf()

ggplot() +
  geom_sf(data = st_geometry(na)) +
  geom_sf(data = taiga, fill = 'red3')

if(! file.exists('figures/taiga-test/taiga-map.png')) {
  ggsave('figures/taiga-test/taiga-map.png', width = 10, height = 6.4,
         units = 'in', dpi = 300, bg = 'white')
}

#' `geom_smooth` of a subset of the first year of data
#' there are too many values > 0.2 in winter
if(FALSE) {
  list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
             pattern = '.nc', full.names = TRUE)[seq(1, 365, by = 10)] %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI')
      r %>%
        crop(., st_transform(taiga, crs(.))) %>%
        mask(st_transform(taiga, crs(.))) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    mutate(doy = yday(date)) %>%
    ggplot(aes(doy, NDVI)) +
    geom_point(alpha = 0.1) +
    geom_smooth(formula = y ~ s(x, k = 5, bs = 'cc'), method = 'gam',
                method.args = list(knots = list(doy = c(0.5, 366.5))),
                fullrange = TRUE) +
    xlim(c(0, 366))
}

# import ndvi data ----
if(file.exists('data/taiga-test/taiga-ndvi.rds')) {
  d <- readRDS('data/taiga-test/taiga-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  d <-
    list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI')
      r %>%
        crop(., st_transform(taiga, crs(.))) %>%
        mask(st_transform(taiga, crs(.))) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r)))) %>%
        return()
    }, .progress = TRUE, seed = NULL) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(month_int = month(date),
           ndvi_clean = case_when(
             # if ndvi < 0.2, there's no reason for concern
             ndvi <= 0.2 ~ ndvi,
             # if ndvi is > 0.2 and too far north, set it to NA
             month_int %in% c(1:4, 11, 12) & y > 60 ~ NA_real_,
             month_int %in% 9:10 & y > 70 ~ NA_real_,
             # if ndvi is > 0.2 but not too far north, keep unchanged
             .default = ndvi))
  plan(sequential)
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 4) %>% # nearest res finer than 0.05x0.05
    crop(st_buffer(taiga, 1e4)) %>% # crop to area
    rast() # convert rasterLayer to SpatRaster
  
  plot(elevs)
  plot(taiga, add = TRUE, col = 'transparent', lwd = 2)
  
  d <- mutate(d,
              elev_m = extract(elevs, select(d, x, y))[, 2],
              year = year(date),
              doy = yday(date))
  range(d$elev_m)
  quantile(d$elev_m, c(0.1, 0.01, 0.001))
  
  d <- mutate(d, elev_m = if_else(elev_m < 0, 0, elev_m))
  
  saveRDS(d, 'data/taiga-test/taiga-ndvi.rds')
}

summary(slice_sample(d, n = 1e7))

# plot the first few NDVI rasters (reducing dataset size to speed plot up)
p_d <- filter(slice(d, 1:(nrow(d) / 44 / 4)), date <= date[1] + 20) %>%
  ggplot() +
  facet_wrap(~ date, nrow = 3) +
  geom_sf(data = taiga, fill = 'grey50') +
  geom_raster(aes(x, y, fill = ndvi_clean)) +
  geom_sf(data = taiga, fill = 'transparent') +
  scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
  ylab(NULL) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
p_d

# plot a sample of the data to check the seasonal trend
layout(t(1:2))
# with original ndvi data
bam(ndvi ~ s(doy, bs = 'cc'), data = slice_sample(d, n = 1e5),
    method = 'REML', knots = list(doy = c(0.5, 366.5))) %>%
  plot(scheme = 2, xlim = c(0, 366), n = 400, resid = TRUE)

# using the cleaned ndvi data
bam(ndvi_clean ~ s(doy, bs = 'cc'), data = slice_sample(d, n = 1e5),
    method = 'REML', knots = list(doy = c(0.5, 366.5))) %>%
  plot(scheme = 2, xlim = c(0, 366), n = 400, resid = TRUE)
layout(1)

# the markov random field is excruciatigly slow for the region. the model
# took over 23 hours to set the bases for a single day of data, so we are
# abboandoning the mrf smooth. the tests below use a single day of data.
if(FALSE) {
  # add cell ID and make list of neighbor cells for all coordinates ----
  # unique coordinates inside taiga
  locs <- d %>%
    select(x, y) %>%
    group_by(x, y) %>%
    slice(1) %>%
    st_as_sf(coords = c('x', 'y')) %>%
    st_set_crs('EPSG:4326') %>%
    filter(., st_as_sf(., coords = c('x', 'y')) %>%
             st_set_crs('EPSG:4326') %>%
             st_intersects(st_transform(taiga, 'EPSG:4326'),
                           sparse = TRUE) %>%
             map_lgl(\(x) length(x) > 0))
  nrow(locs) # all rasters use same coords
  
  p_locs <-
    ggplot() +
    geom_sf(data = taiga) +
    geom_sf(data = locs, alpha = 0.1)
  p_locs
  
  # make a raster of all locations
  r_0 <- locs %>%
    st_coordinates() %>%
    data.frame() %>%
    mutate(z = 0) %>%
    rast()
  
  p_locs +
    geom_raster(aes(x, y), as.data.frame(r_0, xy = TRUE), fill = '#FF000030') +
    labs(x = NULL, y = NULL)
  
  if(file.exists('data/taiga-test/nbs.rds')) {
    nbs <- readRDS('data/taiga-test/nbs.rds')
  } else {
    nbs <- nbs_from_rast(r_0)
    saveRDS(nbs, 'data/taiga-test/nbs.rds')
  }
  
  d <-
    mutate(d,
           cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
             factor(levels = names(nbs)))
  
  if(FALSE) { # check cell IDs
    all.equal(sort(names(nbs)),
              sort(as.character(unique(unlist(nbs)))))
    
    # there is one pixel in the north that has no neighbors
    nbs[which(map_int(nbs, length) == 1 & map_lgl(nbs, \(.n) .n[1] == 0))]
    
    filter(d, cell_id == '501')
  }
  
  # all cell names match the neighbor list names
  # there cannot be any list elements with names that are not in the dataset
  all.equal(sort(as.character(cells(r_0))), sort(names(nbs)))
  
  # all values in neighbor list are in the factor levels
  # this is not crucial, as the model will predict for neighbors with data
  all.equal(sort(levels(d$cell_id)),
            sort(as.character(unique(unlist(nbs[-which(names(nbs) == '501')])))))
  
  all.equal(sort(as.character(unique(d$cell_id))),
            sort(names(nbs)))
  
  # filter to a single day and drop the cell with no neighbors
  z <- filter(d, cell_id != '501') %>%
    filter(date == date[1]) %>%
    filter(! is.na(ndvi_clean)) %>%
    mutate(cell_id = factor(cell_id))
  nbs_z <- nbs[levels(z$cell_id)]
  
  all.equal(sort(names(nbs_z)), sort(levels(z$cell_id)))
  all.equal(sort(levels(z$cell_id)), sort(as.character(unique(unlist(nbs_z)))))
  all.equal(sort(as.character(unique(z$cell_id))), sort(names(nbs_z)))
  all.equal(sort(as.character(unique(z$cell_id))),
            sort(as.character(unique(unlist(nbs_z)))))
  
  # too slow! setting up bases since about 2025-06-19 11:30, and it was
  # still not done at 7:30 the next day (20 hours later)
  system.time(m_test <- bam(
    ndvi_clean ~
      s(cell_id, bs = 'mrf', k = 1000,
        xt = list(nb = nbs_z)),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(nthreads = 10, trace = TRUE)))
  
  # fits in < 25 seconds on the EME linux
  system.time(m_test_ds <- bam(
    ndvi_clean ~
      s(x, y, bs = 'ds', k = 1000),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(nthreads = 10, trace = TRUE)))
  
  # also fits in < 25 seconds on the EME linux
  system.time(m_test_sos <- bam(
    ndvi_clean ~
      s(x, y, bs = 'sos', k = 1000),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(nthreads = 10, trace = TRUE)))
}

# fit a spatially explicit test model with a gaussian family ----
if(file.exists('models/taiga-test/gaussian-gam-ds.rds')) {
  m_gaus_ds <- readRDS('models/taiga-test/gaussian-gam-ds.rds')
} else {
  # fits in ~30 minutes
  m_gaus_ds <- bam(
    ndvi_clean ~
      s(x, y, bs = 'ds', k = 2000) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0.5, 366.5)),
    data = d,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(nthreads = 10, trace = TRUE))
  
  saveRDS(m_gaus_ds, 'models/taiga-test/gaussian-gam-ds.rds')
  p_ds <- draw(m_gaus_ds, rug = FALSE, dist = 0.07)
  ggsave('figures/taiga-test/taiga-ndvi-gaussian-ds-terms.png', p_ds,
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
}

# deviance explained and complexity of the spatial terms are similar
summary(m_gaus_ds)

# # testing data aggregation ----
#' #' *need to add code to clean ndvi data*
# s_res <- 4 # spatial resolution
# t_res <- 4 # temporal resolution
# 
# # import and aggregate the data ----
# if(file.exists('data/taiga-test/taiga-ndvi-t-4-s-4-aggr.rds')) {
#   d_aggr <- readRDS('data/taiga-test/taiga-ndvi-t-4-s-4-aggr.rds')
# } else {
#   if(.Platform$OS.type != 'unix')
#     stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
#   future::availableCores(logical = FALSE)
#   plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
#   d_aggr <-
#     list.files(path = 'data/avhrr-viirs-ndvi/raster-files/', #'/home/shared/NOAA_Files/',
#                pattern = '.nc',
#                full.names = TRUE) %>%
#     future_map(\(fn) {
#       r <- rast(fn, lyr = 'NDVI') # to extract time below
#       r %>%
#         crop(., st_transform(taiga, crs(.))) %>%
#         terra::aggregate(s_res, na.rm = TRUE) %>%
#         mask(st_transform(taiga, crs(.))) %>% # mask after aggregating
#         as.data.frame(xy = TRUE, na.rm = TRUE) %>%
#         mutate(date = as.Date(unique(time(r))))
#     }, .progress = TRUE) %>%
#     bind_rows() %>%
#     as_tibble() %>%
#     rename(ndvi = NDVI) %>%
#     # aggregate temporally
#     mutate(julian = julian(date),
#            central_date = as.Date(julian - (julian %% t_res) + 2)) %>%
#     group_by(central_date, x, y) %>%
#     summarize(doy = yday(central_date),
#               year = year(central_date),
#               ndvi = mean(ndvi, na.rm = TRUE)) %>%
#     ungroup()
#   plan(sequential)
#   
#   # plot the first few NDVI rasters
#   p_d_aggr <-
#     filter(d_aggr, central_date %in% unique(as.character(central_date))[1:21]) %>%
#     ggplot() +
#     facet_wrap(~ as.character(central_date), nrow = 3) +
#     geom_raster(aes(x, y, fill = ndvi)) +
#     geom_sf(data = taiga, fill = 'transparent') +
#     scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
#     ylab(NULL) +
#     scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
#   p_d_aggr
#   
#   # add altitude (get_elev_points results in a json error)
#   elevs_aggr <- d %>%
#     filter(date == first(date)) %>%
#     select(x, y) %>%
#     mutate(z = 1) %>%
#     rast(crs = 'EPSG:4326') %>%
#     aggregate(4, na.rm = TRUE) %>%
#     get_elev_raster(z = 2) %>% # nearest finer res than 0.20x0.20
#     crop(st_buffer(taiga, 1e4))
#   
#   plot(elevs_aggr)
#   plot(taiga, add = TRUE, col = 'transparent')
#   
#   d_aggr <- mutate(d_aggr,
#                    elev_m = extract(elevs_aggr, select(d_aggr, x, y)),
#                    year = year(central_date),
#                    doy = yday(central_date))
#   
#   range(d_aggr$elev_m)
#   quantile(d_aggr$elev_m, c(0.1, 0.01, 0.001))
#   
#   saveRDS(d_aggr, paste0('data/taiga-test/taiga-ndvi-t-',
#                          t_res, '-s-', s_res, '-aggr.rds'))
# }
# 
# # ndvi still has high max values
# summary(d_aggr)
# 
# # plot a comparison of the first 21 rasters
# plot_grid(
#   get_legend(p_d +
#                theme(legend.position = 'top', legend.key.width = rel(2))),
#   plot_grid(p_d + theme(legend.position = 'none'),
#             p_d_aggr + theme(legend.position = 'none'),
#             labels = 'AUTO'),
#   rel_heights = c(1, 15), ncol = 1)
# 
# ggsave('figures/taiga-test/example-rasters-aggregation.png',
#        width = 20, height = 8.5, units = 'in', dpi = 600, bg = 'white')
# 
# # fit a model to show the changes in predictions ----
# # add cell ID and make list of neighbor cells for all coordinates ----
# # unique coordinates inside taiga
# locs_aggr <- d_aggr %>%
#   select(x, y) %>%
#   group_by(x, y) %>%
#   slice(1) %>%
#   st_as_sf(coords = c('x', 'y')) %>%
#   st_set_crs('EPSG:4326') %>%
#   filter(., st_as_sf(., coords = c('x', 'y')) %>%
#            st_set_crs('EPSG:4326') %>%
#            st_transform(crs(taiga)) %>%
#            st_intersects(taiga, sparse = TRUE) %>%
#            map_lgl(\(x) length(x) > 0))
# nrow(locs_aggr) # all rasters use same coords
# 
# p_locs_aggr <-
#   ggplot() +
#   geom_sf(data = taiga) +
#   geom_sf(data = locs_aggr)
# p_locs_aggr
# 
# # make a raster of all locations
# r_0_aggr <- locs_aggr %>%
#   st_coordinates() %>%
#   data.frame() %>%
#   mutate(z = 0) %>%
#   rast()
# 
# p_locs_aggr +
#   geom_raster(aes(x, y), as.data.frame(r_0_aggr, xy = TRUE),
#               fill = '#FF000030') +
#   labs(x = NULL, y = NULL)
# 
# nbs_aggr <- nbs_from_rast(r_0_aggr)
# 
# d_aggr <-
#   mutate(d_aggr,
#          cell_id = cellFromXY(r_0_aggr, xy = as.matrix(tibble(x, y))) %>%
#            factor(levels = names(nbs_aggr))) %>%
#   filter(! is.na(cell_id))
# mean(is.na(d_aggr$cell_id)) # some cell centers fall outside of the polygon
# 
# n_distinct(d_aggr$cell_id)
# length(names(nbs_aggr))
# 
# # all cell names match the neighbor list names
# # there cannot be any list elements with names that are not in the dataset
# all.equal(sort(as.character(cells(r_0_aggr))), sort(names(nbs_aggr)))
# 
# # all values in neighbor list are in the factor levels
# # this is not crucial, as the model will predict for neighbors with data
# all.equal(sort(levels(d_aggr$cell_id)),
#           sort(as.character(unique(unlist(nbs_aggr)))))
# 
# if(file.exists('models/taiga-test/gaussian-gam-ds.rds')) {
#   m_gaus_mrf_aggr <- readRDS('models/taiga-test/gaussian-gam-mrf-aggr.rds')
# } else {
#   m_gaus_mrf_aggr <- bam(
#     ndvi ~
#       s(cell_id, bs = 'mrf', k = 14, xt = list(nb = nbs_aggr)) +
#       s(elev_m, bs = 'cr', k = 5) +
#       s(year, bs = 'cr', k = 10) +
#       s(doy, bs = 'cc', k = 10),
#     family = gaussian(),
#     knots = list(doy = c(0, 1)),
#     data = d_aggr,
#     method = 'fREML',
#     discrete = TRUE)
#   plot_mrf(.model = m_gaus_mrf_aggr, .newdata = d_aggr, .full_model = TRUE)
#   ggsave('figures/taiga-test/taiga-ndvi-gaussian-mrf-aggr-terms.png',
#          width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
#   saveRDS(m_gaus_mrf_aggr, 'models/taiga-test/gaussian-gam-mrf-aggr.rds')
# }
# 
# # make a figure comparing ds gam, mrf gam, and mrf_aggr gam ----
# elevs <- d %>%
#   filter(date == first(date)) %>%
#   select(x, y) %>%
#   mutate(z = 1) %>%
#   rast(crs = 'EPSG:4326') %>%
#   get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
#   crop(st_buffer(taiga, 1e4)) # crop to area near taiga
# 
# elevs_aggr <- d %>%
#   filter(date == first(date)) %>%
#   select(x, y) %>%
#   mutate(z = 1) %>%
#   rast(crs = 'EPSG:4326') %>%
#   aggregate(4, na.rm = TRUE) %>%
#   get_elev_raster(z = 2) %>% # nearest finer res than 0.20x0.20
#   crop(st_buffer(taiga, 1e4))
# 
# gratia::smooths(m_gaus_ds)
# gratia::smooths(m_gaus_mrf)
# 
# # get model predictions ----
# get_preds <- function(nd, space = TRUE) {
#   if(space) {
#     preds <- nd %>%
#       mutate(.,
#              ds_mu =
#                predict(object = m_gaus_ds, newdata = .,
#                        type = 'response', se.fit = FALSE,
#                        terms = c('(Intercept)', 's(x,y)', 's(elev_m)')),
#              mrf_mu =
#                rename(., cell_id = cell_id_fine) %>%
#                predict(object = m_gaus_mrf, newdata = .,
#                        type = 'response', se.fit = FALSE,
#                        terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')),
#              mrfa_mu =
#                rename(., cell_id = cell_id_aggr) %>%
#                predict(object = m_gaus_mrf_aggr, newdata = .,
#                        type = 'response', se.fit = FALSE,
#                        terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')))
#   } else {
#     preds <- nd %>% # new data
#       mutate(
#         ds_mu = predict(object = m_gaus_ds, newdata = .,
#                         type = 'response', se.fit = FALSE,
#                         terms = c('(Intercept)', 's(doy)')),
#         mrf_mu = predict(object = m_gaus_mrf, newdata = .,
#                          type = 'response', se.fit = FALSE,
#                          terms = c('(Intercept)', 's(doy)')),
#         mrfa_mu = predict(object = m_gaus_mrf_aggr, newdata =.,
#                           type = 'response', se.fit = FALSE,
#                           terms = c('(Intercept)', 's(doy)')))
#   }
#   
#   preds <- preds %>%
#     mutate(diff_mu = mrfa_mu - mrf_mu) %>%
#     select(x, y, elev_m, doy, ds_mu:diff_mu)
#   
#   if(space) {
#     # get temporally static maps of estimated variance
#     get_s2 <- function(.model) {
#       .m <- if(.model == 'Duchon spline') {
#         .m <- m_gaus_ds
#       } else if(.model == 'MRF') {
#         .m <- m_gaus_mrf
#       } else if(.model == 'MRF (aggregated)') {
#         .m <- m_gaus_mrf_aggr
#       }
#       
#       # need to select the dataset based on the model but keep (x,y) coords
#       .d <- if(.model == 'Duchon spline' | .model == 'MRF') {
#         .d <- d
#       } else if(.model == 'MRF (aggregated)') {
#         .d <- na.omit(d_aggr)
#       }
#       
#       .d %>%
#         transmute(x, y, e2 = resid(.m)^2) %>%
#         group_by(x, y) %>%
#         summarise(s2 = mean(e2), .groups = 'drop') %>%
#         rast() %>%
#         `crs<-`('EPSG:4326') %>%
#         project(elevs, res = res(elevs)) %>%
#         extract(., select(as.data.frame(elevs, xy = TRUE), 1:2)) %>%
#         select(! ID) %>%
#         bind_cols(select(as.data.frame(elevs, xy = TRUE), 1:2), .) %>%
#         mutate(model = .model) %>%
#         filter(! is.na(s2)) %>%
#         as_tibble() %>%
#         mutate(model = .model) %>%
#         return()
#     }
#     
#     s2 <- bind_rows(map(c('Duchon spline', 'MRF', 'MRF (aggregated)'), get_s2)) %>%
#       pivot_wider(names_from = model, values_from = s2) %>%
#       mutate(diff = `MRF (aggregated)` - `MRF`) %>%
#       pivot_longer(`Duchon spline`:diff, names_to = 'model', values_to = 'value')
#     
#     preds <-
#       bind_rows(pivot_longer(preds, ds_mu:diff_mu,
#                              names_to = c('model', 'param'),
#                              values_to = 'value', names_sep = '_'),
#                 mutate(s2, param = 's2'))
#   } else {
#     preds <- 
#       left_join(
#         preds,
#         tibble(
#           doy = unique(d$doy),
#           ds_s2 = d %>%
#             mutate(e2 = resid(m_gaus_ds)^2) %>%
#             group_by(doy) %>%
#             summarize(s2 = mean(e2, na.rm = TRUE)) %>%
#             pull(s2),
#           mrf_s2 = d %>%
#             mutate(e2 = resid(m_gaus_mrf)^2) %>%
#             group_by(doy) %>%
#             summarize(s2 = mean(e2, na.rm = TRUE)) %>%
#             pull(s2)),
#         by = 'doy') %>%
#       left_join(
#         tibble(doy = unique(d_aggr$doy),
#                mrfa_s2 = d_aggr %>%
#                  na.omit() %>%
#                  mutate(e2 = resid(m_gaus_mrf_aggr)^2) %>%
#                  group_by(doy) %>%
#                  summarize(s2 = mean(e2, na.rm = TRUE)) %>%
#                  pull(s2)),
#         by = 'doy') %>%
#       mutate(diff_s2 = mrfa_s2 - mrf_s2) %>%
#       select(doy, ds_mu:diff_s2) %>%
#       pivot_longer(ds_mu:diff_s2,
#                    names_to = c('model', 'param'), values_to = 'value',
#                    names_sep = '_')
#   }
#   preds %>%
#     mutate(model = case_when(model == 'ds' ~ 'Duchon spline',
#                              model == 'mrf' ~ 'MRF',
#                              model == 'mrfa' ~ 'MRF (aggregated)',
#                              model == 'diff' ~ 'diff',
#                              TRUE ~ model)) %>%
#     return()
# }
# 
# # spatial predictions
# preds_comp_s <-
#   elevs %>%
#   mask(taiga) %>%
#   as.data.frame(xy = TRUE) %>%
#   rename(elev_m = 3) %>%
#   filter(! is.na(elev_m)) %>%
#   mutate(elev_m = if_else(elev_m < 0, 0, elev_m)) %>%
#   mutate(cell_id_fine = factor(cells(r_0, vect(tibble(x, y),
#                                                geom = c('x', 'y')))[, 2],
#                                levels = levels(d$cell_id)),
#          cell_id_aggr = factor(cells(r_0_aggr,
#                                      vect(tibble(x, y),
#                                           geom = c('x', 'y')))[, 2],
#                                levels = levels(d_aggr$cell_id)),
#          year = 0, doy = 0) %>%
#   get_preds(space = TRUE) # add model predictions
# 
# # temporal doy predictions
# preds_comp_t <- tibble(doy = 1:366,
#                        cell_id = intersect(d$cell_id, d_aggr$cell_id)[1],
#                        x = 0, y = 0, elev_m = 0, year = 0) %>%
#   get_preds(space = FALSE)
# 
# # make the final figure ----
# p_comp <-
#   plot_grid(
#     ncol = 2, labels = 'AUTO', rel_widths = c(3, 1.5),
#     rel_heights = c(3, 3, 2, 2, 2),
#     # row 1: map of mean NDVI
#     ggplot(filter(preds_comp_s, param == 'mu', model != 'diff')) +
#       facet_grid(. ~ model) +
#       geom_raster(aes(x, y, fill = value)) +
#       geom_sf(data = taiga, fill = 'transparent', color = 'black') +
#       scale_fill_viridis_c('NDVI', option = 'A') +
#       labs(x = NULL, y = NULL),
#     ggplot(filter(preds_comp_s, param == 'mu', model == 'diff')) +
#       geom_raster(aes(x, y, fill = value)) +
#       geom_sf(data = taiga, fill = 'transparent', color = 'black') +
#       scale_fill_distiller(
#         expression(atop(bold('Difference in'), bold('mean NDVI'))),
#         type = 'div', palette = 5,
#         limits = max(abs(filter(preds_comp_s, model == 'diff',
#                                 param == 'mu')$value)) * c(-1, 1)) +
#       labs(x = NULL, y = NULL),
#     # row 2: map of variance in NDVI
#     filter(preds_comp_s, param == 's2', model != 'diff') %>%
#       group_by(model) %>%
#       mutate(value = value / max(value, na.rm = TRUE)) %>%
#       ggplot() +
#       facet_grid(. ~ model) +
#       geom_raster(aes(x, y, fill = value)) +
#       geom_sf(data = taiga, fill = 'transparent', color = 'black') +
#       scale_fill_viridis_c(expression(bold(s^'2'))) +
#       labs(x = NULL, y = NULL),
#     ggplot(filter(preds_comp_s, param == 's2', model == 'diff')) +
#       geom_raster(aes(x, y, fill = value)) +
#       geom_sf(data = taiga, fill = 'transparent', color = 'black') +
#       scale_fill_distiller(
#         expression(bold(atop('Difference in', 'estimated s'^'2'))),
#         type = 'div', palette = 4) +
#       labs(x = NULL, y = NULL),
#     # row 3: mean NDVI over day of year
#     ggplot(filter(preds_comp_t, param == 'mu', model != 'diff')) +
#       facet_grid(. ~ model) +
#       geom_line(aes(doy, value)) +
#       labs(x = 'Day of year', y = 'Mean NDVI'),
#     ggplot(filter(preds_comp_t, param == 'mu', model == 'diff')) +
#       geom_line(aes(doy, value)) +
#       geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
#       labs(x = 'Day of year', y = 'Difference in mean NDVI'),
#     # row 4: variance in NDVI over day of year
#     ggplot(filter(preds_comp_t, param == 's2', model != 'diff'),
#            aes(doy, value)) +
#       facet_grid(. ~ model) +
#       geom_point(alpha = 0.3) +
#       geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
#                   method.args = list(knots = list(x = c(0.5, 366.5)))) +
#       geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
#       labs(x = 'Day of year',
#            y = expression(bold(paste('Daily mean e'^'2')))),
#     ggplot(filter(preds_comp_t, param == 's2', model == 'diff'),
#            aes(doy, value)) +
#       geom_point(alpha = 0.3) +
#       geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
#                   method.args = list(knots = list(x = c(0.5, 366.5)))) +
#       geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
#       geom_rug(aes(x = x), tibble(x = unique(d$doy)), alpha = 0.1,
#                inherit.aes = FALSE) +
#       labs(x = 'Day of year',
#            y = expression(bold(paste('Difference in daily mean e'^'2')))),
#     # row 5: NDVI over elevation
#     ggplot(filter(preds_comp_s, param == 'mu', model != 'diff') %>%
#              filter(! is.na(value))) +
#       facet_grid(. ~ model) +
#       geom_point(aes(elev_m, value), alpha = 0.1) +
#       geom_rug(aes(x = elev_m), alpha = 0.1) +
#       labs(x = 'Elevation (m)', y = 'Mean NDVI'),
#     ggplot(filter(preds_comp_s, param == 'mu', model == 'diff') %>%
#              filter(! is.na(value))) +
#       geom_point(aes(elev_m, value), alpha = 0.1) +
#       geom_smooth(aes(elev_m, value), formula = y ~ s(x, k = 5),
#                   method = 'gam') +
#       geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
#       geom_rug(aes(x = elev_m), alpha = 0.01) +
#       labs(x = 'Elevation (m)', y = 'Difference in mean NDVI'))
# 
# ggsave('figures/taiga-test/model-comparisons.png',
#        width = 12.5, height = 17, units = 'in', dpi = 300, bg = 'white')
