# specs of my personal laptop:
# 64.0 GB RAM, 13th Gen Intel Core i7-1370P processor, 14 cores
# specs of the EME Linux:
# 2.2 TB RAM, 2 Intel Xeon Platinum 8462Y+ processors, 64 cores
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
source('functions/betals.r') # custom beta location-scale family
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('functions/qr.default_with_LAPACK.R') # for fitting large betals GAM
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')
source('functions/nbs_from_rast.R')
source('functions/is_flagged.R') # to flag bad values based on QA layer
source('functions/decode_qa.R')
assignInNamespace(x = 'qr.default',   # replace `base::qr.default()`
                  value = qr.default, # with the local version sourced above
                  ns = 'base')        # specify the base package

# shapefile of WWF terrestrial ecoregions
sardinia <- read_sf('data/wwf-ecoregions/wwf_terr_ecos.shp') %>%
  st_make_valid() %>% #' necessary for `st_intersection()`
  st_intersection(
    # the bounding box for sardinia (large enough to include all coast)
    tibble(x = c(7.6, 10.4), y = c(38.8, 41.3)) %>%
      st_as_sf(coords = c('x', 'y')) %>%
      st_bbox() %>%
      st_as_sfc() %>%
      st_as_sf() %>%
      st_set_crs('EPSG:4326')) %>%
  st_cast('POLYGON', warn = FALSE) %>%
  st_geometry() %>%
  st_as_sf()

# check QA flags
if(FALSE) {
  r_0 <- list.files(path = 'data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1',
                    pattern = '.nc', full.names = TRUE)[2] %>%
    rast() %>%
    crop(st_transform(sardinia, crs(.)), mask = TRUE)
  
  decoded_0 <- decode_qa(as.numeric(values(r_0$QA)), sensor = 'AVHRR')
  
  # make rasters for each QA parameter
  for(cn in colnames(decoded_0)[3:ncol(decoded_0)]) {
    r_0[[cn]] <- r_0$QA
    values(r_0[[cn]]) <- decoded_0[[cn]]
    r_0[[cn]] <- mask(r_0[[cn]], r_0$QA)
  }
  plot(r_0)
  
  r_0 <- list.files(path = 'data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001/',
                    pattern = '.nc', full.names = TRUE) %>%
    `[`(., length(.)) %>% # take the last raster
    rast() %>%
    crop(st_transform(sardinia, crs(.)), mask = TRUE)
  
  decoded_0 <- decode_qa(as.numeric(values(r_0$QA)), sensor = 'VIIRS')
  
  # make rasters for each QA parameter
  for(cn in colnames(decoded_0)[3:ncol(decoded_0)]) {
    r_0[[cn]] <- r_0$QA
    values(r_0[[cn]]) <- decoded_0[[cn]]
    r_0[[cn]] <- mask(r_0[[cn]], r_0$QA)
  }
  plot(r_0)
  
  rm(r_0, decoded_0)
}

# import ndvi data ----
file_names <- list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                         pattern = '.nc', full.names = TRUE,
                         recursive = TRUE)

if(file.exists('data/sardinia-test/sardinia-ndvi.rds')) {
  d <- readRDS('data/sardinia-test/sardinia-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  d <- file_names %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = c('NDVI', 'QA')) %>%
        crop(., st_transform(sardinia, crs(.))) %>%
        mask(., st_transform(sardinia, crs(.)))
      
      r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI) # drop cloudy pixels
      
      # AVHRR v6 rasters don't have dates, so we need to add them manually
      if(grepl('AVHRR', fn)) {
        .date <- gsub('.*\\N07_AVH13C1.A', '', fn) %>%
          gsub('\\..*', '', .) %>%
          as.Date(format = '%Y%j')
      } else {
        .date <- as.Date(unique(time(r)))
      }
      
      as.data.frame(r$NDVI, xy = TRUE, na.rm = TRUE) %>%
        mutate(date = .date) %>%
        return()
    }, .progress = TRUE, .options = furrr_options(seed = NULL)) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(ndvi_scaled = ndvi_to_01(ndvi))
  plan(sequential)
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    slice(1:1e4) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
    crop(st_buffer(sardinia, 1e4)) # crop to area near sardinia
  
  plot(elevs)
  plot(sardinia, add = TRUE, col = 'transparent')
  
  d <- mutate(d, 
              elev_m = terra::extract(elevs, select(d, x, y)),
              year = year(date),
              doy = yday(date))
  range(d$elev_m)
  quantile(d$elev_m, c(0.1, 0.01, 0.001))
  
  saveRDS(d, 'data/sardinia-test/sardinia-ndvi.rds')
}

# plot the first few NDVI rasters
p_d <- filter(d, date %in% unique(date)[1:21]) %>%
  ggplot() +
  facet_wrap(~ date, nrow = 3) +
  geom_raster(aes(x, y, fill = ndvi)) +
  geom_sf(data = sardinia, fill = 'transparent') +
  scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
  ylab(NULL) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))

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
           st_transform(crs(sardinia)) %>%
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

nbs <- nbs_from_rast(r_0)

d <-
  mutate(d,
         cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs)))

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0))), sort(names(nbs)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
# the values differ because one cell has no neighbors, so it returns 0
# and is not included in the list elements
all.equal(sort(levels(d$cell_id)),
          sort(as.character(unique(unlist(nbs)))))

# test spatial terms with a single day of data ----
if(FALSE) {
  # data for only the first day
  d_0 <- filter(d, date == max(date))
  
  # not all points have data
  ggplot(d_0) +
    geom_raster(aes(x, y, fill = ndvi)) +
    geom_sf(data = locs, color = 'darkorange')
  
  m_mrf_0 <-
    bam(ndvi ~ s(cell_id, bs = 'mrf', k = 200,
                 # need to subset the neighbors to the factor levels
                 xt = list(nb = nbs[levels(d_0$cell_id)])),
        family = gaussian(),
        data = d_0,
        method = 'fREML',
        discrete = TRUE,
        drop.unused.levels = FALSE)
  
  m_ds_0 <- bam(ndvi ~ s(x, y, bs = 'ds', k = 200),
                family = gaussian(),
                data = d_0,
                method = 'fREML',
                discrete = TRUE,
                drop.unused.levels = TRUE)
  
  ggplot() +
    coord_equal() +
    geom_point(aes(predict(m_mrf_0, newdata = d_0),
                   predict(m_ds_0, newdata = d_0)),
               alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    labs(x = 'Fitted values from the cell-level MRF GAM',
         y = 'Fitted values from Duchon-spline GAM')
  
  ggsave(paste0('figures/sardinia-test/duchon-vs-cell-mrf-',
                unique(d_0$date), '.png'),
         width = 5, height = 4, units = 'in', dpi = 600, bg = 'white')
  
  # MRF model is more flexible, gives a better fit, and has good shrinkage
  plot_grid(
    plot_mrf(.model = m_mrf_0, .terms = c('(Intercept)', 's(cell_id)'),
             .newdata = d_0) +
      ggtitle('Cell-level MRF'),
    plot_mrf(.model = m_mrf_0, .terms = c('(Intercept)', 's(cell_id)'),
             .newdata = d_0, .limits = c(NA, NA),
             .pal = viridis::viridis(10)) +
      ggtitle(''),
    d_0 %>%
      mutate(., mu_hat = predict(m_mrf_0, newdata = .)) %>%
      ggplot(aes(ndvi, mu_hat)) +
      coord_equal() +
      geom_point(alpha = 0.2) +
      geom_smooth(method = 'gam', formula = y ~ s(x)) +
      geom_abline(intercept = 0, slope = 1, color = 'red') +
      labs(x = 'Fitted', y = 'Observed'),
    
    draw(m_ds_0, dist = 0.02) & coord_sf(crs = 'EPSG:4326'),
    draw(m_ds_0, dist = 0.02, fun = \(x) x + coef(m_ds_0)['(Intercept)']) &
      coord_sf(crs = 'EPSG:4326') &
      # to specify limits manually
      scale_fill_viridis_c('NDVI', limits = range(fitted(m_mrf_0))),
    d_0 %>%
      mutate(mu_hat = predict(m_ds_0)) %>%
      ggplot(aes(ndvi, mu_hat)) +
      coord_equal() +
      geom_point(alpha = 0.2) +
      geom_smooth(method = 'gam', formula = y ~ s(x)) +
      geom_abline(intercept = 0, slope = 1, color = 'red') +
      labs(x = 'Fitted', y = 'Observed'),
    nrow = 2, rel_widths = c(1, 1, 1.5))
  ggsave(paste0('figures/sardinia-test/duchon-vs-cell-mrf-',
                unique(d_0$date), '-predictions.png'),
         width = 12.5, height = 10, units = 'in', dpi = 300, bg = 'white')
}

# fit a spatially explicit test model with a gaussian family ----
# MRF model is faster and gives more flexible predictions while keeping
# reasonable values
if(all(file.exists(c('models/sardinia-test/gaussian-gam-ds.rds',
                     'models/sardinia-test/gaussian-gam-mrf.rds')))) {
  m_gaus_ds <- readRDS('models/sardinia-test/gaussian-gam-ds.rds')
  m_gaus_mrf <- readRDS('models/sardinia-test/gaussian-gam-mrf.rds')
} else {
  # fits in ~ 15 s
  system.time(
    m_gaus_ds <- bam(
      ndvi ~
        s(x, y, bs = 'ds', k = 200) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = gaussian(),
      knots = list(doy = c(0.5, 366.5)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      nthreads = 10,
      control = gam.control(trace = TRUE))
  )
  
  # fits in ~6 s
  system.time(
    m_gaus_mrf <- bam(
      ndvi ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = gaussian(),
      knots = list(doy = c(0.5, 366.5)),
      data = d,
      method = 'fREML',
      discrete = TRUE,
      control = gam.control(nthreads = 1, trace = TRUE))
  )
  
  # deviance explained and complexity of the spatial terms are similar
  summary(m_gaus_ds)
  summary(m_gaus_mrf)
  
  plot_grid(
    draw(m_gaus_ds, select = 1, rug = FALSE, dist = 0.07) &
      geom_sf(data = sardinia, inherit.aes = FALSE, fill = 'transparent',
              linewidth = 1),
    draw(m_gaus_ds, select = 2, rug = FALSE, dist = 0.07),
    draw(m_gaus_ds, select = 3, rug = FALSE, dist = 0.07),
    draw(m_gaus_ds, select = 4, rug = FALSE, dist = 0.07))
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-ds-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  plot_mrf(.model = m_gaus_mrf, .newdata = d, .full_model = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-mrf-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  saveRDS(m_gaus_ds, 'models/sardinia-test/gaussian-gam-ds.rds')
  saveRDS(m_gaus_mrf, 'models/sardinia-test/gaussian-gam-mrf.rds')
}

# fit a spatially explicit test model with a beta family ----
# fits in ~14 minutes
if(file.exists('models/sardinia-test/beta-gam-mrf.rds')) {
  m_beta_mrf <- readRDS('models/sardinia-test/beta-gam-mrf.rds')
} else {
  system.time(
    m_beta_mrf <- bam(
      ndvi_scaled ~
        s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
        s(elev_m, bs = 'cr', k = 5) +
        s(year, bs = 'cr', k = 10) +
        s(doy, bs = 'cc', k = 10),
      family = betar(),
      knots = list(doy = c(0.5, 366.5)),
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
#' requires `LAPACK = TRUE` in `base::qr.default()`, which can only be done
#' with `assignInNamespace()` because `fixInNamespace()` can't edit the
#' function.
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
      knots = list(doy = c(0.5, 366.5)),
      data = d,
      method = 'REML',
      samfrac = 0.001,
      control = gam.control(nthreads = 10, trace = TRUE))
  )
}

# plots of predicted means and squared residuals for day 1 ----
preds_1 <- d %>%
  filter(date == max(date) - 5) %>%
  mutate(mu_gaus = predict(m_gaus_mrf, newdata = .),
         mu_beta = ndvi_to_11(predict(m_beta_mrf, newdata = .,
                                      type = 'response')))

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
  pivot_longer(cols = c(mu_gaus, mu_beta), names_to = 'model',
               names_prefix = 'mu_', values_to = 'mu') %>%
  mutate(model = case_when(model == 'beta' ~ 'Beta MRF GAM',
                           model == 'gaus' ~ 'Gaussian MRF GAM'),
         e2 = (ndvi - mu)^2)

cowplot::plot_grid(
  # plot estimated mean
  ggplot() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = mu), preds_1_long) +
    geom_sf(data = sardinia, fill = 'transparent') +
    scale_fill_gradientn(expression(paste(hat('\U1D53C'), '(\U1D6CE)')),
                         colors = ndvi_pal, limits = c(-1, 1)) +
    labs(title = paste('Predictions for', d$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_long) +
    geom_sf(data = sardinia, fill = 'transparent') +
    scale_fill_viridis_c(expression(paste(hat('\U1D53C'), '(\U1D6CE - ',
                                          hat('\U1D53C'), '(\U1D6CE)', '\U00B2)')),
      limits = c(0, NA)) +
    labs(title = paste('Squared residuals', d$date[1]), x = NULL, y = NULL),
  labels = 'AUTO', ncol = 1)

ggsave('sardinia-ndvi-mrf-model-agreement-gaussian-beta-predictions.png',
       path = 'figures/sardinia-test', height = 12, width = 8, bg = 'white')

# check agreement with other models
p_fits_mrf <-
  mutate(d,
         gaus = predict(m_gaus_mrf, newdata = d, type = 'response'),
         beta = predict(m_beta_mrf, newdata = d, type = 'response') %>%
           ndvi_to_11()) %>%
  filter(! is.na(cell_id)) %>% #' drop 4% with `NA` `cell_id`
  ggplot() +
  geom_point(aes(x, x), data = tibble(x = 0:1), color = 'transparent') +
  geom_hex(aes(gaus, beta), bins = 75) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Gaussian model (with MRF smooth)',
       y = 'NDVI predicted by Beta model (with MRF smooth)') +
  scale_fill_lapaz(name = 'Count', reverse = TRUE)

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
      knots = list(doy = c(0.5, 366.5)),
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
s_res <- 2 # factor for aggregating spatial resolution
t_res <- 2 # factor for aggregating temporal resolution

# import and aggregate the data ----
if(file.exists(paste0('data/sardinia-test/sardinia-ndvi-t-', t_res,
                      '-s-', s_res, '-aggr.rds'))) {
  d_aggr <- readRDS(paste0('data/sardinia-test/sardinia-ndvi-t-', t_res,
                           '-s-', s_res, '-aggr.rds'))
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  d_aggr <- file_names %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = c('NDVI', 'QA')) %>%
        crop(., st_transform(sardinia, crs(.)), mask = TRUE)
      
      r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI) # drop cloudy pixels
      
      # AVHRR v6 rasters don't have dates, so we need to add them manually
      if(grepl('AVHRR', fn)) {
        .date <- gsub('.*\\N.._AVH13C1.A', '', fn) %>%
          gsub('\\..*', '', .) %>%
          as.Date(format = '%Y%j')
      } else {
        .date <- as.Date(unique(time(r)))
      }
      
      r$NDVI %>%
        terra::aggregate(s_res, na.rm = TRUE) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = .date)
    }, .progress = TRUE, .options = furrr_options(seed = NULL)) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    # aggregate temporally
    mutate(julian = julian(date),
           central_date = as.Date(julian - (julian %% t_res) + t_res / 2)) %>%
    group_by(central_date, x, y) %>%
    summarize(doy = unique(yday(central_date)), # sometimes returnes 2 values
              year = unique(year(central_date)), # sometimes returnes 2 values
              ndvi = mean(ndvi, na.rm = TRUE)) %>%
    ungroup()
  plan(sequential)
  
  # add altitude (get_elev_points results in a json error)
  elevs_aggr <- d %>%
    slice(1:1e4) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 3) %>% # nearest finer res than 0.10x0.10
    crop(st_buffer(sardinia, 1e4))
  
  plot(elevs_aggr)
  plot(sardinia, add = TRUE, col = 'transparent')
  
  d_aggr <- d_aggr %>%
    mutate(elev_m = terra::extract(elevs_aggr, select(d_aggr, x, y)),
           year = year(central_date),
           doy = yday(central_date))
  
  range(d_aggr$elev_m)
  quantile(d_aggr$elev_m, c(0.1, 0.01, 0.001))
  
  saveRDS(d_aggr, paste0('data/sardinia-test/sardinia-ndvi-t-',
                         t_res, '-s-', s_res, '-aggr.rds'))
}

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

# plot a comparison of the first 21 rasters
plot_grid(
  get_legend(p_d +
               theme(legend.position = 'top', legend.key.width = rel(2))),
  p_d + theme(legend.position = 'none'),
            p_d_aggr + theme(legend.position = 'none'),
            labels = c('', 'A', 'B'),
  rel_heights = c(1, 15, 15), ncol = 1)

ggsave('figures/sardinia-test/example-rasters-aggregation.png',
       width = 10, height = 16.5, units = 'in', dpi = 600, bg = 'white')

summary(d_aggr)

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
           st_transform(st_crs(sardinia)) %>%
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

nbs_aggr <- nbs_from_rast(r_0_aggr)

d_aggr <-
  mutate(d_aggr,
         cell_id = cellFromXY(r_0_aggr, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs_aggr)))

mean(is.na(d_aggr$cell_id)) # some cell centers fall outside of the polygon
d_aggr <- filter(d_aggr, ! is.na(cell_id))
n_distinct(d_aggr$cell_id)
length(names(nbs_aggr))

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0_aggr))), sort(names(nbs_aggr)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
# the values differ because one cell has no neighbors, so it returns 0
# and is not included in the list elements
all.equal(sort(levels(d_aggr$cell_id)),
          sort(as.character(unique(unlist(nbs_aggr)))))

if(file.exists('models/sardinia-test/gaussian-gam-ds.rds')) {
  m_gaus_mrf_aggr <- readRDS('models/sardinia-test/gaussian-gam-mrf-aggr.rds')
} else {
  m_gaus_mrf_aggr <- bam(
    ndvi ~
      s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs_aggr)) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0.5, 366.5)),
    data = filter(d_aggr, ! is.na(cell_id)),
    method = 'fREML',
    discrete = TRUE)
  plot_mrf(.model = m_gaus_mrf_aggr, .newdata = d_aggr, .full_model = TRUE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-mrf-aggr-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  saveRDS(m_gaus_mrf_aggr, 'models/sardinia-test/gaussian-gam-mrf-aggr.rds')
}

# make a figure comparing ds gam, mrf gam, and mrf_aggr gam ----
elevs <- d %>%
  slice(1:1e4) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
  crop(st_buffer(sardinia, 1e4)) # crop to area near sardinia

elevs_aggr <- d %>%
  slice(1:1e4) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 3) %>% # nearest finer res than 0.10x0.10
  crop(st_buffer(sardinia, 1e4))

gratia::smooths(m_gaus_ds)
gratia::smooths(m_gaus_mrf)

# get model predictions ----
get_preds <- function(nd, space = TRUE) {
  if(space) {
    preds <- nd %>%
      mutate(.,
             ds_mu =
               predict(object = m_gaus_ds, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(x,y)', 's(elev_m)')),
             mrf_mu =
               rename(., cell_id = cell_id_fine) %>%
               predict(object = m_gaus_mrf, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')),
             mrfa_mu =
               rename(., cell_id = cell_id_aggr) %>%
               predict(object = m_gaus_mrf_aggr, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')))
  } else {
    preds <- nd %>% # new data
      mutate(
        ds_mu = predict(object = m_gaus_ds, newdata = .,
                        type = 'response', se.fit = FALSE,
                        terms = c('(Intercept)', 's(doy)')),
        mrf_mu = predict(object = m_gaus_mrf, newdata = .,
                         type = 'response', se.fit = FALSE,
                         terms = c('(Intercept)', 's(doy)')),
        mrfa_mu = predict(object = m_gaus_mrf_aggr, newdata = .,
                          type = 'response', se.fit = FALSE,
                          terms = c('(Intercept)', 's(doy)')))
  }
  
  preds <- preds %>%
    mutate(diff_mu = mrfa_mu - mrf_mu) %>%
    select(x, y, elev_m, doy, ds_mu:diff_mu)
  
  if(space) {
    # get temporally static maps of estimated variance
    get_s2 <- function(.model) {
      .m <- if(.model == 'Duchon spline') {
        .m <- m_gaus_ds
      } else if(.model == 'MRF') {
        .m <- m_gaus_mrf
      } else if(.model == 'MRF (aggregated)') {
        .m <- m_gaus_mrf_aggr
      }
      
      # need to select the dataset based on the model but keep (x,y) coords
      .d <- if(.model == 'Duchon spline' | .model == 'MRF') {
        .d <- d
      } else if(.model == 'MRF (aggregated)') {
        .d <- na.omit(d_aggr)
      }
      
      .d %>%
        transmute(x, y, e2 = (predict(.m, newdata = .) - ndvi)^2) %>%
        group_by(x, y) %>%
        summarise(s2 = mean(e2), .groups = 'drop') %>%
        rast() %>%
        `crs<-`('EPSG:4326') %>%
        project(elevs, res = res(elevs)) %>%
        mask(st_transform(sardinia, crs(elevs)), touches = FALSE) %>%
        terra::extract(., select(as.data.frame(elevs, xy = TRUE), 1:2)) %>%
        select(! ID) %>%
        bind_cols(select(as.data.frame(elevs, xy = TRUE), 1:2), .) %>%
        mutate(model = .model) %>%
        filter(! is.na(s2)) %>%
        as_tibble() %>%
        mutate(model = .model) %>%
        return()
    }
    
    s2 <- bind_rows(map(c('Duchon spline', 'MRF', 'MRF (aggregated)'), get_s2)) %>%
      pivot_wider(names_from = model, values_from = s2) %>%
      mutate(diff = `MRF (aggregated)` - `MRF`) %>%
      pivot_longer(`Duchon spline`:diff, names_to = 'model', values_to = 'value')
    
    preds <-
      bind_rows(pivot_longer(preds, ds_mu:diff_mu,
                             names_to = c('model', 'param'),
                             values_to = 'value', names_sep = '_'),
                mutate(s2, param = 's2'))
  } else {
    preds <- preds %>%
      mutate(
        ds_s2 = d %>%
          mutate(e2 = (predict(m_gaus_ds, newdata = .) - ndvi)^2) %>%
          group_by(doy) %>%
          summarize(s2 = mean(e2)) %>%
          pull(s2),
        mrf_s2 = d %>%
          mutate(e2 = (predict(m_gaus_mrf, newdata = .) - ndvi)^2) %>%
          na.omit() %>%
          group_by(doy) %>%
          summarize(s2 = mean(e2)) %>%
          pull(s2),
        mrfa_s2 = d_aggr %>%
          na.omit() %>%
          mutate(e2 = (predict(m_gaus_mrf_aggr, newdata = .) - ndvi)^2) %>%
          group_by(doy) %>%
          summarize(s2 = mean(e2)) %>%
          pull(s2),
        diff_s2 = mrfa_s2 - mrf_s2) %>%
      select(doy, ds_mu:diff_s2) %>%
      pivot_longer(ds_mu:diff_s2,
                   names_to = c('model', 'param'), values_to = 'value',
                   names_sep = '_')
  }
  preds %>%
    mutate(model = case_when(model == 'ds' ~ 'Duchon spline',
                             model == 'mrf' ~ 'MRF',
                             model == 'mrfa' ~ 'MRF (aggregated)',
                             model == 'diff' ~ 'diff',
                             TRUE ~ model)) %>%
    return()
}

# spatial predictions
preds_comp_s <-
  elevs %>%
  mask(sardinia) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(elev_m = 3) %>%
  filter(! is.na(elev_m)) %>%
  mutate(cell_id_fine = factor(cells(r_0, vect(tibble(x, y),
                                               geom = c('x', 'y')))[, 2],
                               levels = levels(d$cell_id)),
         cell_id_aggr = factor(cells(r_0_aggr,
                                     vect(tibble(x, y),
                                          geom = c('x', 'y')))[, 2],
                               levels = levels(d_aggr$cell_id)),
         year = 0, doy = 0, elev_m = if_else(elev_m > 0, elev_m, 0)) %>%
  get_preds(space = TRUE) # add model predictions

# temporal doy predictions
preds_comp_t <- tibble(doy = 1:366, x = 0, y = 0, elev_m = 0, year = 0,
                       cell_id = factor(intersect(d$cell_id,
                                                  d_aggr$cell_id)[1])) %>%
  get_preds(space = FALSE)

# make the final figure ----
p_comp <-
  plot_grid(
    ncol = 2, labels = 'AUTO', rel_widths = c(3, 1.5),
    rel_heights = c(3, 3, 2, 2, 2),
    # row 1: map of mean NDVI
    ggplot(filter(preds_comp_s, param == 'mu', model != 'diff')) +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = sardinia, fill = 'transparent', color = 'black') +
      scale_fill_gradientn('Mean NDVI', colors = ndvi_pal, limits = c(-1, 1)) +
      labs(x = NULL, y = NULL),
    ggplot(filter(preds_comp_s, param == 'mu', model == 'diff')) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = sardinia, fill = 'transparent', color = 'black') +
      scale_fill_distiller(
        expression(atop(bold('Difference in'), bold('mean NDVI'))),
        type = 'div', palette = 5,
        limits = max(abs(filter(preds_comp_s, model == 'diff',
                                param == 'mu')$value), na.rm = TRUE) *
          c(-1, 1)) +
      labs(x = NULL, y = NULL),
    # row 2: map of variance in NDVI
    filter(preds_comp_s, param == 's2', model != 'diff') %>%
      group_by(model) %>%
      mutate(value = if_else(value > 0.05, 0.05, value)) %>%
      ggplot() +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = sardinia, fill = 'transparent', color = 'black') +
      scale_fill_viridis_c(expression(bold('s'^'2'))) +
      labs(x = NULL, y = NULL),
    ggplot(filter(preds_comp_s, param == 's2', model == 'diff')) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = sardinia, fill = 'transparent', color = 'black') +
      scale_fill_distiller(
        expression(bold(atop('Difference', 'in s'^'2'))),
        type = 'div', palette = 4,
        limits = max(abs(filter(preds_comp_s, model == 'diff',
                                param == 's2')$value), na.rm = TRUE) *
          c(-1, 1)) +
      labs(x = NULL, y = NULL),
    # row 3: mean NDVI over day of year
    ggplot(filter(preds_comp_t, param == 'mu', model != 'diff')) +
      facet_grid(. ~ model) +
      geom_line(aes(doy, value)) +
      labs(x = 'Day of year', y = 'Mean NDVI'),
    ggplot(filter(preds_comp_t, param == 'mu', model == 'diff')) +
      geom_line(aes(doy, value)) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      labs(x = 'Day of year', y = 'Difference in mean NDVI'),
    # row 4: variance in NDVI over day of year
    ggplot(filter(preds_comp_t, param == 's2', model != 'diff'),
           aes(doy, value)) +
      facet_grid(. ~ model) +
      geom_point(alpha = 0.3) +
      geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
                  method.args = list(knots = list(x = c(0.5, 366.5)))) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      labs(x = 'Day of year',
           y = expression(bold(paste('Daily mean e'^'2')))),
    ggplot(filter(preds_comp_t, param == 's2', model == 'diff'),
           aes(doy, value)) +
      geom_point(alpha = 0.3) +
      geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
                  method.args = list(knots = list(x = c(0.5, 366.5)))) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      geom_rug(aes(x = unique(d$doy)), alpha = 0.1, inherit.aes = FALSE) +
      labs(x = 'Day of year',
           y = expression(bold(paste('Difference in daily mean e'^'2')))),
    # row 5: NDVI over elevation
    ggplot(filter(preds_comp_s, param == 'mu', model != 'diff') %>%
             filter(! is.na(value))) +
      facet_grid(. ~ model) +
      geom_point(aes(elev_m, value), alpha = 0.1) +
      geom_smooth(aes(elev_m, value), formula = y ~ s(x, k = 5),
                  method = 'gam') +
      geom_rug(aes(x = elev_m), alpha = 0.1) +
      labs(x = 'Elevation (m)', y = 'Mean NDVI'),
    ggplot(filter(preds_comp_s, param == 'mu', model == 'diff') %>%
             filter(! is.na(value))) +
      geom_point(aes(elev_m, value), alpha = 0.1) +
      geom_smooth(aes(elev_m, value), formula = y ~ s(x, k = 5),
                  method = 'gam') +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      geom_rug(aes(x = elev_m), alpha = 0.01) +
      labs(x = 'Elevation (m)', y = 'Difference in mean NDVI'))

ggsave('figures/sardinia-test/model-comparisons.png', p_comp,
       width = 12.5, height = 20, units = 'in', dpi = 300, bg = 'white')

# plots not used
if(FALSE) {
  filter(preds_comp_s, param == 'mu', model != 'diff') %>%
    pivot_wider(values_from = value, names_from = model) %>%
    filter(! is.na(`MRF (aggregated)`), ! is.na(`MRF`)) %>%
    ggplot() +
    geom_point(aes(MRF, `MRF (aggregated)`), alpha = 0.1) +
    geom_abline(intercept = 0, slope = 1, color = 'red')
  
  filter(preds_comp_s, param == 'mu', model != 'diff') %>%
    pivot_wider(values_from = value, names_from = model) %>%
    filter(! is.na(`MRF (aggregated)`), ! is.na(`MRF`)) %>%
    ggplot() +
    geom_hex(aes(MRF, `MRF (aggregated)`)) +
    geom_abline(intercept = 0, slope = 1, color = 'red') +
    labs(x = 'MRF', y = 'Aggregated MRF') +
    scale_fill_lapaz(name = 'Count', reverse = TRUE)
}

# make a figure comparing the effects of aggregation on prediction ----
pred_vs_pred <- tibble(
  mrf = d_aggr %>%
    mutate(.,
           cell_id = extract(slice(d, 1, .by = c(x, y)) %>%
                               select(x, y, cell_id) %>%
                               rast(),select(d_aggr, x, y))[[2]]) %>%
    predict(m_gaus_mrf, newdata = ., type= 'response'),
  mrf_aggr = fitted(m_gaus_mrf_aggr)) %>%
  filter(! is.na(mrf))

ggplot(pred_vs_pred, aes(mrf, mrf_aggr)) +
  geom_hex(bins = 75) +
  geom_smooth(method = 'gam', formula = y ~ s(x, k = 5),
              color = 'cornflowerblue') +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'MRF model with full dataset',
       y = 'MRF model with aggregated dataset') +
  scale_fill_lapaz(name = 'Count', reverse = TRUE) +
  theme(legend.position = 'inside', legend.position.inside = c(0.1, 0.85))

ggsave('figures/sardinia-test/pred-vs-pred.png', width = 8, height = 6,
       units = 'in', dpi = 300, bg = 'white')
