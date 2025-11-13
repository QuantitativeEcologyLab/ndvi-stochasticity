# specs of the EME Linux:
# 2.2 TB RAM, 2 Intel Xeon Platinum 8462Y+ processors, 64 cores
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
source('functions/ndvi-palette.R')
source('functions/decode_qa.R') # for cleaning rasters
source('functions/bit_to_int.R') # for cleaning rasters
source('functions/is_flagged.R') # for cleaning rasters
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R') # get_legend() from cowplot v. 1.1.3 fails
source('functions/nbs_from_rast.R') # gives a list of neighboring cells
source('analysis/figures/000-robinson-objects.R') # for Robinson projection

# pick a northern polygon
taiga <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(ecoregion == 'Northwest Territories taiga' |
                  ecoregion == 'Muskwa-Slave Lake forests') %>%
  st_union() %>%
  st_as_sf() %>%
  # drop islands that have long >> 0
  st_transform(crs(rast('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc')))

world <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_union() %>%
  st_transform(robinson_crs)

ggplot() +
  geom_sf(data = bounds, fill = 'white') +
  geom_sf(data = world) +
  geom_sf(data = st_transform(taiga, robinson_crs), fill = 'red',
          color = 'red4') +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_void()

ggsave('figures/taiga-test/taiga-map.png', width = 8, height = 4,
       units = 'in', dpi = 600, bg = 'white')

# some very large NDVI values at the edge of the preliminary cleaning ----
#' `geom_smooth` of a subset of the first year of data
#' there are too many NDVI values > 0 in winter
decode_qa(bit_to_int('0000000000000000'), sensor = 'AVHRR')$cloudy
decode_qa(bit_to_int('0100000000000000'), sensor = 'AVHRR')$cloudy
decode_qa(bit_to_int('0000000000000000'), sensor = 'AVHRR')$cloud_shadow
decode_qa(bit_to_int('0010000000000000'), sensor = 'AVHRR')$cloud_shadow

decode_qa(bit_to_int('0000000000000000'), sensor = 'VIIRS')$cloud_state
decode_qa(bit_to_int('1000000000000000'), sensor = 'VIIRS')$cloud_state
decode_qa(bit_to_int('0100000000000000'), sensor = 'VIIRS')$cloud_state
decode_qa(bit_to_int('1100000000000000'), sensor = 'VIIRS')$cloud_state
decode_qa(bit_to_int('0000000000000000'), sensor = 'VIIRS')$cloud_shadow
decode_qa(bit_to_int('0010000000000000'), sensor = 'VIIRS')$cloud_shadow

# some quick tests for why it is worth dropping cloudy pixels
if(FALSE) {
  list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1',
             pattern = 'A198221.', full.names = TRUE) %>%
    map(\(x) {
      .r <- rast(x, lyr = c('NDVI', 'QA')) %>%
        crop(taiga, mask = TRUE)
      .r$QA <- is_flagged(.r$QA, 1) # check for cloud cover (bits 0 & 1)
      return(.r)
    }) %>%
    rast() %>%
    plot()
  
  z <-
    list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/',
               pattern = '.nc', full.names = TRUE)[seq(1, 365, by = 10)] %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = c('NDVI', 'QA'))
      # AVHRR v6 rasters don't have dates, so we need to add them manually
      if(grepl('AVHRR', fn)) {
        .date <- gsub('.*\\N.._AVH13C1.A', '', fn) %>%
          gsub('\\..*', '', .) %>%
          as.Date(format = '%Y%j')
      } else {
        .date <- as.Date(unique(time(r)))
      }
      
      r %>%
        crop(., st_transform(taiga, crs(.))) %>%
        mask(st_transform(taiga, crs(.))) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = .date) %>%
        return()
    }, .progress = TRUE, .options = furrr_options(seed = NULL)) %>%
    bind_rows() %>%
    slice_sample(n = 1e5) %>%
    mutate(doy = yday(date),
           month = month(date),
           # find confidently or probably cloudy pixels (pos starts from 0)
           cloudy = map_lgl(QA, \(.qa) is_flagged(.qa, flag_position = 1)),
           # find pixels with cloud shadow
           shadow = map_lgl(QA, \(.qa) is_flagged(.qa, 2))) %>%
    filter(! (y > 60 & NDVI > 0.2 & month %in% c(1:4, 11, 12)))
  
  m_z_0 <- bam(NDVI ~ s(doy, bs = 'cc', k = 10), data = z,
               method = 'fREML', discrete = TRUE,
               knots = list(doy = c(0.5, 366.5)))
  m_z_1 <- bam(NDVI ~ s(doy, bs = 'cc', k = 10), data = filter(z, !cloudy),
               method = 'fREML', discrete = TRUE,
               knots = list(doy = c(0.5, 366.5)))
  m_z_2 <- bam(NDVI ~ s(doy, bs = 'cc', k = 10), data = filter(z, !shadow),
               method = 'fREML', discrete = TRUE,
               knots = list(doy = c(0.5, 366.5)))
  
  pred_ndvi <- function(m) {
    tibble(doy = 1:366) %>%
      mutate(., ndvi = predict(m, newdata = .))
  }
  
  plot(ndvi ~ doy, pred_ndvi(m_z_0), type = 'l', ylim = c(-0.1, 1))
  # removing cloudy pixels increases NDVI but dramatically reduces dataset
  lines(ndvi ~ doy, pred_ndvi(m_z_1), col = 2)
  # removing cloud shadows barely has an effect
  lines(ndvi ~ doy, pred_ndvi(m_z_2), col = 3)
  
  rm(z, m_z_0, m_z_1, m_z_2, pred_ndvi)
  
  # showing a few more examples
  z <- list.files(path = 'data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/',
                  pattern = 'A1982217', full.names = TRUE) %>%
    rast(lyr = c('NDVI', 'QA')) %>%
    crop(taiga, mask = TRUE) %>%
    as.data.frame(xy = TRUE) %>%
    filter(! is.na(NDVI)) %>%
    bind_cols(., decode_qa(.$QA, sensor = 'AVHRR')) %>%
    as_tibble() %>%
    mutate(across(cloudy:polar_flag, factor))
  
  summary(z)
  summarize(z, prop = n(), .by = cloudy) %>%
    mutate(prop = prop / sum(prop))
  
  #' `cloud_shadow` and `snow_ice` are not helpful
  plot_grid(
    ggplot(z, aes(x, y, fill = NDVI)) +
      coord_sf(crs = 'EPSG:4326') +
      geom_raster() +
      scale_fill_viridis_c(),
    # scale_fill_gradientn(colours = ndvi_pal, limits = c(-1, 1)),
    ggplot(z, aes(x, y, fill = cloudy)) +
      coord_sf(crs = 'EPSG:4326') +
      geom_raster() +
      scale_fill_bright(),
    ggplot(z, aes(x, y, fill = cloud_shadow)) +
      coord_sf(crs = 'EPSG:4326') +
      geom_raster() +
      scale_fill_bright(),
    ggplot(z, aes(x, y, fill = polar_flag)) +
      coord_sf(crs = 'EPSG:4326') +
      geom_raster() +
      scale_fill_bright())
  
  # dropping cloudy pixels drops some good values but mostly bad ones
  ggplot(z, aes(NDVI)) +
    facet_wrap(~ cloudy) +
    geom_histogram()
  
  rm(z)
}

# import ndvi data ----
file_names <- list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                         pattern = '.nc', full.names = TRUE,
                         recursive = TRUE)

if(file.exists('data/taiga-test/taiga-ndvi.rds')) {
  d <- readRDS('data/taiga-test/taiga-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  
  # takes < 30 minutes on EME linux
  d <-
    future_map(file_names, \(fn) {
      r <- rast(fn, lyr = c('NDVI', 'QA')) %>%
        crop(., st_transform(taiga, crs(.)), mask = TRUE)
      
      r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI) # drop cloudy pixels
      
      # AVHRR v6 rasters don't have dates, so we need to add them manually
      if(grepl('AVHRR', fn)) {
        .date <- gsub('.*\\N.._AVH13C1.A', '', fn) %>%
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
    mutate(month_int = month(date),
           # drop unrealistically high NDVI values in cold periods
           #' see `analysis/005-check-high-lat-data.R` for more info
           ndvi = case_when(
             # if ndvi <= 0, there's no reason for concern
             ndvi <= 0 ~ ndvi,
             # if ndvi is > 0:
             #    and too far north, set it to NA
             month_int == 1 & y > 60 ~ NA_real_,
             month_int == 2 & y > 65 ~ NA_real_,
             month_int == 3 & y > 70 ~ NA_real_,
             month_int == 9 & y > 80 ~ NA_real_,
             month_int == 10 & y > 69 ~ NA_real_,
             month_int == 11 & y > 62 ~ NA_real_,
             month_int == 12 & y > 59 ~ NA_real_,
             #    and not too far north, keep unchanged
             .default = ndvi))
  plan(sequential)
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    slice(1:1e4) %>%
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
  
  saveRDS(d, 'data/taiga-test/taiga-ndvi.rds')
}

summary(slice_sample(d, n = 1e7))

# filtering by lat and NDVI dropped some extreme NDVI values (~2%)
mean(is.na(d$ndvi))

# use cleaned NDVI data and drop NAs
d <- filter(d, ! is.na(ndvi))

# plot the first few NDVI rasters (reducing dataset size to speed plot up)
p_d <- filter(slice(d, 1:(nrow(d) / 44 / 4)), date <= date[1] + 20) %>%
  ggplot() +
  facet_wrap(~ date, nrow = 3) +
  geom_sf(data = taiga, fill = 'grey50') +
  geom_raster(aes(x, y, fill = ndvi)) +
  geom_sf(data = taiga, fill = 'transparent') +
  scale_x_continuous(NULL, breaks = c(-135, -120)) +
  ylab(NULL) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
p_d

# plot a sample of the data to check the seasonal trend
if(FALSE) {
  bam(ndvi ~ s(doy, bs = 'cc'), data = slice_sample(d, n = 1e5),
      method = 'fREML', discrete = TRUE, 
      knots = list(doy = c(0.5, 366.5))) %>%
    plot(scheme = 2, xlim = c(0, 366), n = 400, resid = TRUE)
}

#' the markov random field is excruciatingly slow for the region. the model
#' took over 24 hours to set the bases for a single day of data, so we are
#' abandoning the `mrf` smooth. the tests below use a single day of data.
if(FALSE) {
  # add cell ID and make list of neighbor cells for all coordinates
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
    geom_raster(aes(x, y), as.data.frame(r_0, xy = TRUE),
                fill = '#FF000030') +
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
    filter(! is.na(ndvi)) %>%
    mutate(cell_id = factor(cell_id))
  nbs_z <- nbs[levels(z$cell_id)]
  
  all.equal(sort(names(nbs_z)), sort(levels(z$cell_id)))
  all.equal(sort(levels(z$cell_id)), sort(as.character(unique(unlist(nbs_z)))))
  all.equal(sort(as.character(unique(z$cell_id))), sort(names(nbs_z)))
  all.equal(sort(as.character(unique(z$cell_id))),
            sort(as.character(unique(unlist(nbs_z)))))
  
  # too slow! setting up bases since took over 24 hours (process killed
  # before setup was completed)
  system.time(m_test <- bam(
    ndvi ~
      s(cell_id, bs = 'mrf', k = 1000,
        xt = list(nb = nbs_z)),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(trace = TRUE)))
  
  # fits in < 25 seconds on the EME linux
  system.time(m_test_ds <- bam(
    ndvi ~
      s(x, y, bs = 'ds', k = 1000),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(trace = TRUE)))
  
  # also fits in < 25 seconds on the EME linux
  system.time(m_test_sos <- bam(
    ndvi ~
      s(y, x, bs = 'sos', k = 1000),
    family = gaussian(),
    data = z,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(trace = TRUE)))
}

# fit a spatially explicit test model with a gaussian family ----
if(file.exists('models/taiga-test/gaussian-gam-sos.rds')) {
  m_gaus <- readRDS('models/taiga-test/gaussian-gam-sos.rds')
} else {
  # fits in < 15 minutes
  m_gaus <- bam(
    ndvi ~
      s(y, x, bs = 'sos', k = 200) +
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
  
  saveRDS(m_gaus, 'models/taiga-test/gaussian-gam-sos.rds')
  p_sos <-
    plot_grid(
      draw(m_gaus, select = 1, rug = FALSE, dist = 0.07) &
        geom_sf(data = taiga, inherit.aes = FALSE, fill = 'transparent',
                linewidth = 0.5),
      draw(m_gaus, select = 2, rug = FALSE, dist = 0.07),
      draw(m_gaus, select = 3, rug = FALSE, dist = 0.07),
      draw(m_gaus, select = 4, rug = FALSE, dist = 0.07))
  ggsave('figures/taiga-test/taiga-ndvi-gaussian-sos-terms.png', p_sos,
         width = 9, height = 8, units = 'in', dpi = 300, bg = 'white')
}

# testing data aggregation ----
s_res <- 2 # factor for aggregating spatial resolution
t_res <- 2 # factor for aggregating temporal resolution
AGGR <- paste0('data/taiga-test/taiga-ndvi-t-', t_res, '-s-', s_res,
               '-aggr.rds')

# import and aggregate the data ----
if(file.exists(AGGR)) {
  d_aggr <- readRDS(AGGR)
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  
  # takes ~ 2 hours on EME linux
  d_aggr <-
    future_map(file_names, \(fn) {
      .r <- rast(fn, lyr = c('NDVI', 'QA')) %>%
        crop(., st_transform(taiga, crs(.)), mask = TRUE)
      
      # drop cloudy pixels
      .r$NDVI <- ifel(is_flagged(.r$QA, 1), NA, .r$NDVI)
      
      # AVHRR v6 rasters don't have dates, so we need to add them manually
      if(grepl('AVHRR', fn)) {
        .date <- gsub('.*\\N.._AVH13C1.A', '', fn) %>%
          gsub('\\..*', '', .) %>%
          as.Date(format = '%Y%j')
      } else {
        .date <- as.Date(unique(time(.r)))
      }
      
      terra::aggregate(.r$NDVI, s_res, na.rm = TRUE) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = .date)
    }, .progress = TRUE, .options = furrr_options(seed = NULL)) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(month_int = month(date),
           # drop unrealistically high NDVI values in cold periods
           #' see `analysis/005-check-high-lat-data.R` for more info
           ndvi = case_when(
             # if ndvi <= 0, there's no reason for concern
             ndvi <= 0 ~ ndvi,
             # if ndvi is > 0:
             #    and too far north, set it to NA
             month_int == 1 & y > 60 ~ NA_real_,
             month_int == 2 & y > 65 ~ NA_real_,
             month_int == 3 & y > 70 ~ NA_real_,
             month_int == 9 & y > 80 ~ NA_real_,
             month_int == 10 & y > 69 ~ NA_real_,
             month_int == 11 & y > 62 ~ NA_real_,
             month_int == 12 & y > 59 ~ NA_real_,
             #    and not too far north, keep unchanged
             .default = ndvi),
           julian = julian(date),
           group = julian - (julian %% t_res)) %>%
    # aggregate temporally
    group_by(group, x, y) %>%
    summarize(central_date = as.Date(group + t_res / 2),
              doy = yday(central_date),
              year = year(central_date),
              ndvi = mean(ndvi, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(! is.na(ndvi))
  plan(sequential)
  
  # add altitude (get_elev_points results in a json error)
  elevs_aggr <- d_aggr %>%
    slice(1:1e4) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 3) %>% # nearest finer res than 0.10x0.10
    crop(st_buffer(taiga, 1e4)) %>%
    rast()
  
  plot(elevs_aggr)
  plot(taiga, add = TRUE, col = 'transparent')
  
  d_aggr <- mutate(d_aggr,
                   elev_m = extract(elevs_aggr, select(d_aggr, x, y)))
  
  # very few elevations < 0
  range(d_aggr$elev_m)
  quantile(d_aggr$elev_m, c(0.1, 0.01, 0.001))
  mean(d_aggr$elev_m < 0)
  
  saveRDS(d_aggr, AGGR)
}

# plot the first few NDVI rasters
p_d_aggr <-
  filter(slice(d_aggr, 1:(nrow(d) / 44 / 4)),
         central_date <= min(central_date) + 40) %>%
  ggplot() +
  facet_wrap(~ as.character(central_date), nrow = 3) +
  geom_sf(data = taiga, fill = 'grey50') +
  geom_raster(aes(x, y, fill = ndvi)) +
  geom_sf(data = taiga, fill = 'transparent') +
  scale_x_continuous(NULL, breaks = c(-135, -120)) +
  ylab(NULL) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
p_d_aggr

# plot a comparison of the first 21 rasters
p_d_both <-
  plot_grid(
    get_legend(p_d +
                 theme(legend.position = 'top', legend.key.width = rel(2))),
    p_d + theme(legend.position = 'none'),
              p_d_aggr + theme(legend.position = 'none'),
    labels = c('', 'A', 'B'),
    rel_heights = c(1, 10, 10), ncol = 1)

ggsave('figures/taiga-test/example-rasters-aggregation.png', p_d_both,
       width = 10, height = 12, units = 'in', dpi = 600, bg = 'white')

if(file.exists('models/taiga-test/gaussian-gam-sos-aggr.rds')) {
  m_gaus_aggr <- readRDS('models/taiga-test/gaussian-gam-sos-aggr.rds')
} else {
  # fits in ~ 6 minutes
  m_gaus_aggr <- bam(
    ndvi ~
      s(y, x, bs = 'sos', k = 200) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0.5, 366.5)),
    data = d_aggr,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10,
    control = gam.control(trace = TRUE))
  
  saveRDS(m_gaus_aggr, 'models/taiga-test/gaussian-gam-sos-aggr.rds')
  p_sos_aggr <-
    plot_grid(
      draw(m_gaus_aggr, select = 1, rug = FALSE, dist = 0.07) &
        geom_sf(data = taiga, inherit.aes = FALSE, fill = 'transparent',
                linewidth = 0.5),
      draw(m_gaus_aggr, select = 2, rug = FALSE, dist = 0.07),
      draw(m_gaus_aggr, select = 3, rug = FALSE, dist = 0.07),
      draw(m_gaus_aggr, select = 4, rug = FALSE, dist = 0.07))
  ggsave('figures/taiga-test/taiga-ndvi-gaussian-sos-aggr-terms.png',
         p_sos_aggr, width = 9, height = 8, units = 'in', dpi = 300,
         bg = 'white')
}

# make a figure comparing ds gam, mrf gam, and mrf_aggr gam ----
elevs <- d %>%
  filter(date == first(date)) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
  crop(taiga) # crop to area near taiga

gratia::smooths(m_gaus)
gratia::smooths(m_gaus_aggr)

# get model predictions ----
get_preds <- function(nd, space = TRUE) {
  # get model predictions
  if(space) {
    preds <- nd %>%
      mutate(full_mu =
               predict(object = m_gaus, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(y,x)', 's(elev_m)')),
             aggr_mu =
               predict(object = m_gaus_aggr, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(y,x)', 's(elev_m)')))
  } else {
    preds <- nd %>% # new data
      mutate(
        full_mu = predict(object = m_gaus, newdata = .,
                          type = 'response', se.fit = FALSE,
                          terms = c('(Intercept)', 's(doy)')),
        aggr_mu = predict(object = m_gaus_aggr, newdata = .,
                          type = 'response', se.fit = FALSE,
                          terms = c('(Intercept)', 's(doy)')))
  }
  
  # calculate difference between predictions
  preds <- preds %>%
    mutate(diff_mu = aggr_mu - full_mu) %>%
    select(x, y, elev_m, doy, full_mu:diff_mu)
  
  if(space) {
    # get temporally static maps of estimated variance
    get_s2 <- function(.model) {
      .m <- if(.model == 'full') {
        .m <- m_gaus
      } else if(.model == 'aggr') {
        .m <- m_gaus_aggr
      }
      
      # need to select the dataset based on the model but keep (x,y) coords
      .d <- if(.model == 'full') {
        .d <- d
      } else if(.model == 'aggr') {
        .d <- d_aggr
      }
      
      .d %>%
        transmute(x, y, e2 = resid(.m)^2) %>%
        group_by(x, y) %>%
        summarise(s2 = mean(e2), .groups = 'drop') %>%
        rast() %>%
        `crs<-`('EPSG:4326') %>%
        project(elevs, res = res(elevs)) %>% #' project to same as `elevs`
        # extract the s2 values and add xy coords
        extract(., select(as.data.frame(elevs, xy = TRUE), 1:2)) %>%
        select(! ID) %>%
        bind_cols(select(as.data.frame(elevs, xy = TRUE), 1:2), .) %>%
        mutate(model = .model) %>%
        filter(! is.na(s2)) %>%
        as_tibble() %>%
        return()
    }
    
    s2 <- bind_rows(map(c('full', 'aggr'), get_s2)) %>%
      pivot_wider(names_from = model, values_from = s2) %>%
      mutate(diff = aggr - full) %>%
      pivot_longer(full:diff, names_to = 'model', values_to = 'value')
    
    preds <-
      bind_rows(pivot_longer(preds, full_mu:diff_mu,
                             names_to = c('model', 'param'),
                             values_to = 'value', names_sep = '_'),
                mutate(s2, param = 's2'))
  } else {
    preds <-
      left_join(
        preds,
        tibble(
          doy = unique(d$doy),
          full_s2 = d %>%
            mutate(e2 = resid(m_gaus)^2) %>%
            group_by(doy) %>%
            summarize(s2 = mean(e2, na.rm = TRUE)) %>%
            pull(s2),
          aggr_s2 = d_aggr %>%
            mutate(e2 = resid(m_gaus_aggr)^2) %>%
            group_by(doy) %>%
            summarize(s2 = mean(e2, na.rm = TRUE)) %>%
            pull(s2)),
        by = 'doy') %>%
      mutate(diff_s2 = aggr_s2 - full_s2) %>%
      select(doy, full_mu:diff_s2) %>%
      pivot_longer(full_mu:diff_s2,
                   names_to = c('model', 'param'), values_to = 'value',
                   names_sep = '_')
  }
  preds %>%
    mutate(model = case_when(model == 'full' ~ 'Full dataset',
                             model == 'aggr' ~ 'Aggregated dataset',
                             model == 'diff' ~ 'diff',
                             TRUE ~ model)) %>%
    return()
}

# spatial predictions
preds_comp_s <-
  elevs %>%
  mask(st_transform(taiga, st_crs(elevs))) %>%
  as.data.frame(xy = TRUE, na.rm = TRUE) %>%
  rename(elev_m = 3) %>%
  filter(! is.na(elev_m)) %>%
  mutate(year = 0, doy = 0) %>%
  get_preds(space = TRUE) %>% # add model predictions
  mutate(model = factor(model, levels = c('diff', 'Full dataset',
                                          'Aggregated dataset')))

# temporal doy predictions
preds_comp_t <- tibble(doy = 1:366,
                       x = 0, y = 0, elev_m = 0, year = 0) %>%
  get_preds(space = FALSE) %>% # add model predictions
  mutate(model = factor(model, levels = c('diff', 'Full dataset',
                                          'Aggregated dataset')))

# make the final figure ----
# clamping some values to make trends more visible
preds_comp_s %>%
  filter(model == 'diff') %>%
  ggplot() +
  facet_wrap(~ param, scales = 'free') +
  geom_histogram(aes(value))

preds_comp_s %>%
  filter(model != 'diff', param == 's2') %>%
  ggplot() +
  facet_wrap(~ model, scales = 'free', ncol = 1) +
  geom_histogram(aes(value))

preds_comp_s %>%
  filter(model != 'diff', param == 's2') %>%
  summarize(q99 = quantile(value, 0.99, na.rm = TRUE), .by = model)

preds_comp_s_0 <- preds_comp_s

p_comp <-
  plot_grid(
    ncol = 2, labels = 'AUTO', rel_widths = c(2, 1.25),
    rel_heights = c(3, 3, 2, 2, 2),
    # row 1: map of mean NDVI
    ggplot(filter(preds_comp_s, param == 'mu', model != 'diff')) +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = taiga, fill = 'transparent', color = 'black') +
      scale_fill_viridis_c('NDVI', option = 'A') +
      labs(x = NULL, y = NULL),
    ggplot(filter(preds_comp_s, param == 'mu', model == 'diff')) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = taiga, fill = 'transparent', color = 'black') +
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
      mutate(value = value / quantile(value, 0.975, na.rm = TRUE),
             value = if_else(value > 1, 1, value)) %>%
      ggplot() +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = taiga, fill = 'transparent', color = 'black') +
      scale_fill_viridis_c(expression(bold(s^'2')), limits = c(0, 1)) +
      labs(x = NULL, y = NULL),
    filter(preds_comp_s, param == 's2', model == 'diff') %>%
      mutate(max_value = quantile(abs(value), 0.975, na.rm = TRUE),
             value = if_else(abs(value) > max_value,
                             max_value * sign(value),
                             value)) %>%
      ggplot() +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = taiga, fill = 'transparent', color = 'black') +
      scale_fill_distiller(
        expression(bold(atop('Difference in', 'estimated s'^'2'))),
        type = 'div', palette = 4) +
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
      geom_smooth(formula = y ~ s(x, bs = 'cc', k = 20), method = 'gam',
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
      geom_rug(aes(x = x), tibble(x = unique(d$doy)), alpha = 0.1,
               inherit.aes = FALSE) +
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

ggsave('figures/taiga-test/model-comparisons.png', p_comp,
       width = 12.5, height = 17, units = 'in', dpi = 300, bg = 'white')
