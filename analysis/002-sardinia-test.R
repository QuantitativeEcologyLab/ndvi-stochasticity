# specs of my personal laptop:
# 64.0 GB RAM, 13th Gen Intel Core i7-1370P processor, 14 cores
# specs of the EME Linux:
# 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
library('sf')        # for shapefiles
library('terra')     # for rasters
library('elevatr')   # for digital elevation models
library('dplyr')     # for data wrangling
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('mgcv')      # for GAMs
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids 
library('gratia')    # for fancy plots of GAMs
library('ggplot2')   # for fanct plots
source('functions/betals.r') # custom beta location-scale family
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('functions/qr.default_with_LAPACK.R')
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')
assignInNamespace(x = 'qr.default',   # replace `base::qr.default()`
                  value = qr.default, # with the local version sourced above
                  ns = 'base')        # specify the base package

# shapefile of the world's terrestrial ecoregions
sardinia <- read_sf('data/world-ecosystems/data/commondata/data0/tnc_terr_ecoregions.shp') %>%
  st_make_valid() %>%
  st_intersection(
    # the bounding box for sardinia (large enough to include all coast)
    tibble(x = c(7.6, 10.4), y = c(38.8, 41.3)) %>%
      st_as_sf(coords = c('x', 'y')) %>%
      st_bbox() %>%
      st_as_sfc() %>%
      st_as_sf() %>%
      st_set_crs('EPSG:4326')) %>%
  st_cast('POLYGON') %>%
  st_geometry() %>%
  st_as_sf()

# import ndvi data ----
if(file.exists('data/sardinia-test/sardinia-ndvi.rds')) {
  d <- readRDS('data/sardinia-test/sardinia-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(50, availableCores(logical = FALSE) - 2))
  d <-
    list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI')
      r %>%
        crop(., st_transform(sardinia, crs(.))) %>%
        mask(st_transform(sardinia, crs(.))) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(ndvi_scaled = ndvi_to_01(ndvi))
  plan(sequential)
  
  # plot the first few NDVI rasters
  p_d <- filter(d, date %in% unique(date)[1:21]) %>%
    ggplot() +
    facet_wrap(~ date, nrow = 3) +
    geom_raster(aes(x, y, fill = ndvi)) +
    geom_sf(data = sardinia, fill = 'transparent') +
    scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
    ylab(NULL) +
    scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
    crop(st_buffer(sardinia, 1e4)) # crop to area near sardinia
  
  plot(elevs)
  plot(sardinia, add = TRUE, col = 'transparent')
  
  d <- mutate(d, 
              elev_m = extract(elevs, select(d, x, y)),
              year = year(date),
              doy = yday(date))
  range(d$elev_m)
  quantile(d$elev_m, c(0.1, 0.01, 0.001))
  
  d <- mutate(d, elev_m = if_else(elev_m < 0, 0, elev_m))
  
  saveRDS(d, 'data/sardinia-test/sardinia-ndvi.rds')
}
summary(d)

# add cell ID and make list of neighbor cells for all coordinates ----
# unique coordinates inside sardinia
locs <- d %>%
  select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  filter(., st_as_sf(., coords = c('x', 'y')) %>%
           st_set_crs('EPSG:4326') %>%
           st_intersects(sardinia, sparse = TRUE) %>%
           map_lgl(\(x) length(x) > 0))
nrow(locs) # all rasters use same coords

p_locs <-
  ggplot() +
  geom_sf(data = sardinia) +
  geom_sf(data = locs)
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

nbs <-
  adjacent(r_0, cells = cells(r_0), directions = 8, include = TRUE) %>%
  as.data.frame() %>%
  transmute(ref_cell = V1, # first column is the starting cell
            # add the 8 surrounding neighbors
            adjacent = map(1:n(), \(i) {
              .z <- c(V2[i], V3[i], V4[i], V5[i], V6[i], V7[i], V8[i], V9[i])
              
              .values <- map_lgl(.z, \(.cell_id) {
                if(is.nan(.cell_id)) {
                  return(NA)
                } else {
                  return(r_0[.cell_id]$z[1])
                }})
              
              .z <- .z[which(! is.na(.values))]
              
              if(length(.z) == 0) {
                return(0)
              } else {
                return(as.character(.z))
              }
            }))
names(nbs$adjacent) <- nbs$ref_cell
nbs <- nbs$adjacent

d <-
  mutate(d,
         cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs)))

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0))), sort(names(nbs)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
all.equal(sort(levels(d$cell_id)),
          sort(as.character(unique(unlist(nbs)))))

# test spatial terms with a single day of data ----
# data for only the first day
d_1981_06_25 <- filter(d, date == '1981-06-25')

# not all points have data
ggplot(d_1981_06_25) +
  geom_raster(aes(x, y, fill = ndvi)) +
  geom_sf(data = locs, color = 'darkorange')

m_mrf_1981_06_25 <-
  bam(ndvi ~ s(cell_id, bs = 'mrf', k = 200,
               # need to subset the neighbors to those in the data
               xt = list(nb = nbs[unique(d_1981_06_25$cell_id)])),
      family = gaussian(),
      data = d_1981_06_25,
      method = 'fREML',
      discrete = TRUE,
      drop.unused.levels = TRUE,
      control = gam.control(trace = TRUE))

m_ds_1981_06_25 <- bam(ndvi ~ s(x, y, bs = 'ds'),
                       family = gaussian(),
                       data = d_1981_06_25,
                       method = 'fREML',
                       discrete = TRUE,
                       drop.unused.levels = TRUE,
                       control = gam.control(trace = TRUE))
ggplot() +
  coord_equal() +
  geom_point(aes(fitted(m_mrf_1981_06_25), fitted(m_ds_1981_06_25)),
             alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'Fitted values from the cell-level MRF GAM',
       y = 'Fitted values from Douchon-spline GAM')
ggsave('figures/sardinia-test/douchon-vs-cell-mrf-1981-06-25.png',
       width = 5, height = 4, units = 'in', dpi = 600, bg = 'white')

# MRF model is more flexible, gives a better fit, and has good shrinkage
plot_grid(
  plot_mrf(.model = m_mrf_1981_06_25, .terms = c('(Intercept)', 's(cell_id)'),
           .newdata = d_1981_06_25) +
    ggtitle('Cell-level MRF'),
  plot_mrf(.model = m_mrf_1981_06_25, .terms = c('(Intercept)', 's(cell_id)'),
           .newdata = d_1981_06_25, .limits = c(NA, NA),
           .pal = viridis::viridis(10)) +
    ggtitle(''),
  d_1981_06_25 %>%
    mutate(mu_hat = fitted(m_mrf_1981_06_25)) %>%
    ggplot(aes(ndvi, mu_hat)) +
    coord_equal() +
    geom_point(alpha = 0.2) +
    geom_smooth(method = 'gam', formula = y ~ s(x)) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    labs(x = 'Fitted', y = 'Observed'),
  
  draw(m_ds_1981_06_25, dist = 0.02),
  draw(m_ds_1981_06_25, dist = 0.02, fun = \(x) x + coef(m_ds_1981_06_25)['(Intercept)']) &
    # to specify limits manually
    scale_fill_viridis_c('NDVI', limits = range(fitted(m_mrf_1981_06_25))),
  d_1981_06_25 %>%
    mutate(mu_hat = fitted(m_ds_1981_06_25)) %>%
    ggplot(aes(ndvi, mu_hat)) +
    coord_equal() +
    geom_point(alpha = 0.2) +
    geom_smooth(method = 'gam', formula = y ~ s(x)) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    labs(x = 'Fitted', y = 'Observed'),
  nrow = 2)
ggsave('figures/sardinia-test/douchon-vs-cell-mrf-1981-06-25-predictions.png',
       width = 15, height = 10, units = 'in', dpi = 300, bg = 'white')

# fit a spatially explicit test model with a gaussian family ----
# MRF model is faster and gives more flexible predictions while keeping
# reasonable values
if(all(file.exists(c('models/sardinia-test/gaussian-gam-ds.rds',
                     'models/sardinia-test/gaussian-gam-mrf.rds')))) {
  m_gaus_ds <- readRDS('models/sardinia-test/gaussian-gam-ds.rds')
  m_gaus_mrf <- readRDS('models/sardinia-test/gaussian-gam-mrf.rds')
} else {
  # on personal laptop: fits in 34 s
  # on EME: fits in 30 s
  system.time(
    m_gaus_ds <- bam(
      ndvi ~
        s(x, y, bs = 'ds', k = 200) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = gaussian(),
      knots = list(doy = c(0, 1)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      nthreads = 10,
      control = gam.control(trace = TRUE))
  )
  
  # on personal laptop: fits in 10 s
  # on EME: fits in 30 s
  system.time(
    m_gaus_mrf <- bam(
      ndvi ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = gaussian(),
      knots = list(doy = c(0, 1)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(nthreads = 1, trace = TRUE))
  )
  
  # deviance explained and complexity of the spatial terms are similar
  summary(m_gaus_ds)
  summary(m_gaus_mrf)
  
  draw(m_gaus_ds, rug = FALSE, dist = 0.015)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-ds-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  plot_mrf(.model = m_gaus_mrf, .newdata = d, .full_model = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-mrf-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  saveRDS(m_gaus_ds, 'models/sardinia-test/gaussian-gam-ds.rds')
  saveRDS(m_gaus_mrf, 'models/sardinia-test/gaussian-gam-mrf.rds')
}

# fit a spatially explicit test model with a beta family ----
# on personal laptop: mrf model fits in ~ 21 minutes
#                     ds model fits in ~ 25 miutes
# on EME linux: mrf model fits in ~ 9.0 minutes
#                ds model fits in ~ 8.8 miutes
if(file.exists('models/sardinia-test/beta-gam-mrf.rds')) {
  m_beta_mrf <- readRDS('models/sardinia-test/beta-gam-mrf.rds')
} else {
  system.time(
    m_beta_mrf <- bam(
      ndvi_scaled ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        #s(x, y, bs = 'ds', k = 200) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = betar(),
      knots = list(doy = c(0, 1)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(trace = TRUE))
  )
  saveRDS(m_beta_mrf, 'models/sardinia-test/beta-gam-mrf.rds')
  
  plot_mrf(m_beta_mrf, .newdata = d, .rug = FALSE, .full_model = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  summary(m_beta_mrf)
}

# fit a spatially explicit test model with a betals family ----
#' requires `LAPACK = TRUE` in `base::qr.default()`, which can only be done with
#' `assignInNamespace()` because `fixInNamespace()` can't edit the function
#' fitting time is >> 144 hours = 6 days on EME Linux (3 newtons in 6 days)
#' uses too much RAM (> 68.4 GB) and swap memory (total: > 716 GB)
if(FALSE) {
  system.time(
    m_betals <- gam(list(
      ndvi_scaled ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      ~ s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10)),
      family = betals(),
      knots = list(doy = c(0, 1)),
      data = d,
      method = 'REML',
      samfrac = 0.001,
      control = gam.control(nthreads = 10, trace = TRUE))
  )
  saveRDS(m_betals, 'models/sardinia-test/betals-gam.rds')
  
  draw(m_betals, rug = FALSE)
  ggsave('figures/sardinia-test/sardinia-ndvi-betals-terms.png',
         width = 9, height = 12, units = 'in', dpi = 300, bg = 'white')
  
  summary(m_betals)
}

# plots of predicted means and squared residuals for day 1 ----
preds_1 <- d %>%
  mutate(mu_gaus = fitted(m_gaus_mrf),
         mu_beta = ndvi_to_11(fitted(m_beta_mrf))) %>%
  filter(date == min(date))

ggplot(preds_1, aes(mu_beta, mu_gaus)) +
  geom_point(alpha = 0.1) +
  geom_smooth(method = 'gam', formula = y ~ s(x)) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Beta MRF model (with constant scale parameter)',
       y = 'NDVI predicted by Gaussian MRF model (with constant variance)')
ggsave(paste0('sardinia-ndvi-mrf-model-agreement-gaussian-beta-',
              unique(preds_1$date), '.png'),
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

preds_1_long <- preds_1 %>%
  tidyr::pivot_longer(cols = c(mu_gaus, mu_beta), names_to = 'model',
                      names_prefix = 'mu_', values_to = 'mu') %>%
  mutate(model = case_when(model == 'beta' ~ 'Beta MRF GAM',
                           model == 'gaus' ~ 'Gaussian MRF GAM'),
         e2 = (ndvi - mu)^2)

cowplot::plot_grid(
  # plot estimated mean
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = mu), preds_1_long) +
    scale_fill_viridis_c(expression(hat('\U1D6CE')),
                         option = 'A') +
    labs(title = paste('Predictions for', d$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_long) +
    scale_fill_viridis_c(expression(
      paste('(', hat('\U1D6CE'), ' - \U1D53C(', hat('\U1D6CE'), '\U00B2)')),
      limits = c(0, NA)) +
    labs(title = paste('Squared residuals', d$date[1]), x = NULL, y = NULL),
  labels = 'AUTO', ncol = 1)

ggsave('sardinia-ndvi-mrf-model-agreement-gaussian-beta-predictions.png',
       path = 'figures/sardinia-test', height = 12, width = 8, bg = 'white')

# check agreement with other models
p_fits_mrf <-
  tibble(gaus = predict(m_gaus_mrf, newdata = d, type = 'response'),
         beta = predict(m_beta_mrf, newdata = d, type = 'response') %>%
           ndvi_to_11()) %>%
  ggplot() +
  geom_point(aes(x, x), data = tibble(x = 0:1), color = 'transparent') +
  geom_point(aes(gaus, beta), alpha = 0.01) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Gaussian model (with MRF smooth)',
       y = 'NDVI predicted by Beta model (with MRF smooth)')

ggsave('sardinia-ndvi-model-agreement-mrf.png', plot = p_fits_mrf,
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

# testing complexity of ti(doy, space) ----
#' can fit `ti()` terms with MRF smooths exactly as with other bases
#' *MRF ti fit is flat -- may be because of lack of trends in the data?*
#' on personal laptop: without `ti(doy, mrf)` fits in ~15 s
#' on personal laptop: with `ti(doy, mrf)` fits in ~25 s
if(file.exists('models/sardinia-test/gaus-gam-ti-mrf-gam.rds')) {
  m_gaus_ti_mrf <- readRDS('models/sardinia-test/gaus-gam-ti-mrf-gam.rds')
} else {
  system.time(
    m_gaus_ti_mrf <- bam(
      ndvi ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10) +
        ti(doy, year, bs = c('cc', 'cr'), k = c(5, 5)) +
        ti(doy, cell_id, bs = c('cc', 'mrf'), k = c(5, 10),
           xt = list(nb = nbs)),
      family = gaussian(),
      knots = list(doy = c(0, 1)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(nthreads = 1, trace = TRUE))
  )
  saveRDS(m_gaus_ti_mrf, 'models/sardinia-test/gaus-gam-ti-mrf-gam.rds')
  
  plot_mrf(.model = m_gaus_ti_mrf, .newdata = d, .full_model = TRUE, ti_surfaces = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-mrf-ti-terms.png',
         width = 13.5, height = 8, units = 'in', dpi = 300, bg = 'white')
  plot_mrf(.model = m_gaus_ti_mrf, .newdata = d, .full_model = TRUE, ti_surfaces = FALSE)
}

summary(m_gaus_ti_mrf)

# testing data aggregation ----
s_res <- 4 # spatial resolution
t_res <- 4 # temporal resolution

# import and aggregate the data ----
if(file.exists('data/sardinia-test/sardinia-ndvi-t-4-s-4-aggr.rds')) {
  d_aggr <- readRDS('data/sardinia-test/sardinia-ndvi-t-4-s-4-aggr.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(30, availableCores(logical = FALSE) - 2))
  d_aggr <-
    list.files(path = 'data/avhrr-viirs-ndvi/raster-files/', #'/home/shared/NOAA_Files/',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI') # to extract time below
      r %>%
        crop(., st_transform(sardinia, crs(.))) %>%
        terra::aggregate(s_res, na.rm = TRUE) %>%
        mask(st_transform(sardinia, crs(.))) %>% # mask after aggregating
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    # aggregate temporally
    mutate(julian = julian(date),
           central_date = as.Date(julian - (julian %% t_res) + 2)) %>%
    group_by(central_date, x, y) %>%
    summarize(doy = yday(central_date),
              year = year(central_date),
              ndvi = mean(ndvi, na.rm = TRUE)) %>%
    ungroup()
  plan(sequential)
  
  # plot the first few NDVI rasters
  p_d_aggr <-
    filter(d_aggr, central_date %in% unique(as.character(central_date))[1:21]) %>%
    ggplot() +
    facet_wrap(~ as.character(central_date), nrow = 3) +
    geom_raster(aes(x, y, fill = ndvi)) +
    geom_sf(data = sardinia, fill = 'transparent') +
    scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
    ylab(NULL) +
    scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
  p_d_aggr
  
  # add altitude (get_elev_points results in a json error)
  elevs_aggr <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    aggregate(4, na.rm = TRUE) %>%
    get_elev_raster(z = 2) %>% # nearest finer res than 0.20x0.20
    crop(st_buffer(sardinia, 1e4))
  
  plot(elevs_aggr)
  plot(sardinia, add = TRUE, col = 'transparent')
  
  d_aggr <- mutate(d_aggr, 
                   elev_m = extract(elevs_aggr, select(d_aggr, x, y)),
                   year = year(central_date),
                   doy = yday(central_date))
  
  range(d_aggr$elev_m)
  quantile(d_aggr$elev_m, c(0.1, 0.01, 0.001))
  
  d_aggr <- mutate(d_aggr, elev_m = if_else(elev_m < 0, 0, elev_m))
  
  saveRDS(d_aggr, paste0('data/sardinia-test/sardinia-ndvi-t-',
                         t_res, '-s-', s_res, '-aggr.rds'))
}
summary(d_aggr)

# plot a comparison of the first 21 rasters
plot_grid(
  get_legend(p_d +
               theme(legend.position = 'top', legend.key.width = rel(2))),
  plot_grid(p_d + theme(legend.position = 'none'),
            p_d_aggr + theme(legend.position = 'none'),
            labels = 'AUTO'),
  rel_heights = c(1, 15), ncol = 1)

ggsave('figures/sardinia-test/example-rasters-aggregation.png',
       width = 20, height = 8.5, units = 'in', dpi = 600, bg = 'white')

# fit a model to show the changes in predictions ----
# add cell ID and make list of neighbor cells for all coordinates ----
# unique coordinates inside sardinia
locs_aggr <- d_aggr %>%
  select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  filter(., st_as_sf(., coords = c('x', 'y')) %>%
           st_set_crs('EPSG:4326') %>%
           st_intersects(sardinia, sparse = TRUE) %>%
           map_lgl(\(x) length(x) > 0))
nrow(locs_aggr) # all rasters use same coords

p_locs_aggr <-
  ggplot() +
  geom_sf(data = sardinia) +
  geom_sf(data = locs_aggr)
p_locs_aggr

# make a raster of all locations
r_0_aggr <- locs_aggr %>%
  st_coordinates() %>%
  data.frame() %>%
  mutate(z = 0) %>%
  rast()

p_locs_aggr +
  geom_raster(aes(x, y), as.data.frame(r_0_aggr, xy = TRUE),
              fill = '#FF000030') +
  labs(x = NULL, y = NULL)

nbs_aggr <-
  adjacent(r_0_aggr, cells = cells(r_0_aggr), directions = 8,
           include = TRUE) %>%
  as.data.frame() %>%
  transmute(ref_cell = V1, # first column is the starting cell
            # add the 8 surrounding neighbors
            adjacent = map(1:n(), \(i) {
              .z <- c(V2[i], V3[i], V4[i], V5[i], V6[i], V7[i], V8[i], V9[i])
              
              .values <- map_lgl(.z, \(.cell_id) {
                if(is.nan(.cell_id)) {
                  return(NA)
                } else {
                  return(r_0_aggr[.cell_id]$z[1])
                }})
              
              .z <- .z[which(! is.na(.values))]
              
              if(length(.z) == 0) {
                return(0)
              } else {
                return(as.character(.z))
              }
            }))
names(nbs_aggr$adjacent) <- nbs_aggr$ref_cell
nbs_aggr <- nbs_aggr$adjacent

d_aggr <-
  mutate(d_aggr,
         cell_id = cellFromXY(r_0_aggr, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs_aggr)))
mean(is.na(d_aggr$cell_id)) # some cell centers fall outside of the polygon
n_distinct(d_aggr$cell_id)

length(names(nbs_aggr))
n_distinct(d_aggr$cell_id)

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0_aggr))), sort(names(nbs_aggr)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
all.equal(sort(levels(d_aggr$cell_id)),
          sort(as.character(unique(unlist(nbs_aggr)))))


#' **HERE**
if(file.exists('models/sardinia-test/gaussian-gam-ds.rds')) {
  m_gaus_mrf_aggr <- readRDS('models/sardinia-test/gaussian-gam-mrf-aggr.rds')
} else {
  system.time(
    m_gaus_mrf_aggr <- bam(
      ndvi ~
        s(cell_id, bs = 'mrf', k = 60, xt = list(nb = nbs_aggr)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = gaussian(),
      knots = list(doy = c(0, 1)),
      data = filter(d_aggr, ! is.na(cell_id)),
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(nthreads = 1, trace = TRUE))
  )
  plot_mrf(.model = m_gaus_mrf_aggr, .newdata = d_aggr, .full_model = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-mrf-aggr-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  saveRDS(m_gaus_mrf_aggr, 'models/sardinia-test/gaussian-gam-mrf-aggr.rds')
}

# make a figure comparing ds gam, mrf gam, and mrf_aggr gam ----

#' *include:*
#' row 1: DS
#' row 2: MRF
#' row 3: MRF aggr
#' row 4:
#'    - comparison of preds vs preds (MRF vs MRF aggr)
#'    - raster of difference of spatial terms (use raster for finer model)
#'    - var estimates of two models

# plot_grid(
#   plot_mrf(.model = m_mrf_1981_06_25, .terms = c('(Intercept)', 's(cell_id)'),
#            .newdata = d_1981_06_25) +
#     ggtitle('Cell-level MRF'),
#   plot_mrf(.model = m_mrf_1981_06_25, .terms = c('(Intercept)', 's(cell_id)'),
#            .newdata = d_1981_06_25, .limits = c(NA, NA),
#            .pal = viridis::viridis(10)) +
#     ggtitle(''),
#   d_1981_06_25 %>%
#     mutate(mu_hat = fitted(m_mrf_1981_06_25)) %>%
#     ggplot(aes(ndvi, mu_hat)) +
#     coord_equal() +
#     geom_point(alpha = 0.2) +
#     geom_smooth(method = 'gam', formula = y ~ s(x)) +
#     geom_abline(intercept = 0, slope = 1, color = 'red') +
#     labs(x = 'Fitted', y = 'Observed'),
#   
#   draw(m_ds_1981_06_25, dist = 0.02),
#   draw(m_ds_1981_06_25, dist = 0.02, fun = \(x) x + coef(m_ds_1981_06_25)['(Intercept)']) &
#     # to specify limits manually
#     scale_fill_viridis_c('NDVI', limits = range(fitted(m_mrf_1981_06_25))),
#   d_1981_06_25 %>%
#     mutate(mu_hat = fitted(m_ds_1981_06_25)) %>%
#     ggplot(aes(ndvi, mu_hat)) +
#     coord_equal() +
#     geom_point(alpha = 0.2) +
#     geom_smooth(method = 'gam', formula = y ~ s(x)) +
#     geom_abline(intercept = 0, slope = 1, color = 'red') +
#     labs(x = 'Fitted', y = 'Observed'),
#   nrow = 2)
# 
# ggsave('figures/sardinia-test/model-comparisons.png',
#        width = 15, height = 20, units = 'in', dpi = 300, bg = 'white')
