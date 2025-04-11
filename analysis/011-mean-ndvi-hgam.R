library('dplyr')   # for data wrangling
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('sf')      # for simple fleature objects
library('terra')   # for rasters
library('gratia')  # for plotting GAMs
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/add_nb.R')

ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# resolution is (0.05 degrees)^2 = (3 degree minutes)^2
# create raster of ecoregion
r_eco <- rast('data/ecoregions/wwf-ecoregions.tif')

# create raster of polygon ID
r_poly_id <- rast('data/ecoregions/ecoregion-polygon-id.tif')

# get raster of elevations
r_elev <- rast('data/elev-raster.tif')

# create raster of distance from coast
r_dist <- rast('data/distance-from-coast-m.tif')

# set up the data for the model
if(file.exists('data/hbam-ndvi-data.rds')) {
  d <- readRDS('data/hbam-ndvi-data.rds')
} else {
  d <- readRDS('data/ndvi-global-15-day-average.rds') %>%
    #' use `exact_extract()` for more data?
    mutate(
      doy = lubridate::yday(central_date),
      poly_id = exactextract::exact_extract(r_poly_id, data.frame(x, y))[, 2], # 1 is pixel ID
      wwf_ecoregion = exactextract::exact_extract(r_eco, data.frame(x, y))[, 2],
      # distance_coast_m = exactextract::exact_extract(r_dist, data.frame(x, y))[, 2],
      elevation_m = exactextract::exact_extract(r_elev, data.frame(x, y))[, 2]) %>%
    # convert strings to factors
    mutate(wwf_ecoregion = paste(wwf_ecoregion, if_else(y < 0, 'S', 'N')),
           wwf_ecoregion = factor(wwf_ecoregion))
  
  apply(d, \(.row) any(is.na(.row))) %>%
    mean()
  
  d <- na.omit(d)
  
  # some excessively low elevations for coastal pixels sea level
  # occurs mostly near islands; leaving values as they are because it helps
  # inform about the lower NDVI due to the coast
  if(FALSE) {
    plot(st_geometry(ecoregions))
    d %>%
      filter(elevation_m < -20) %>%
      points(y ~ x, ., col = 'red', pch = 19)
  }
  
  d <- mutate(d, elevation_m = if_else(elevation_m < 0, 0, elevation_m))
  
  saveRDS(d, 'data/hbam-ndvi-data.rds')
}

# create list of neighbors ----
# drop polygons with no data or neighbors
ecoregions$in_data <- ecoregions$poly_id %in% unique(d$poly_id)
filter(ecoregions, ! in_data & n_neigh == 0)

if(FALSE) {
  ggplot(ecoregions, aes(log10(area_km2), in_data)) + geom_point()
  
  ggplot(ecoregions, aes(log10(area_km2), n_neigh)) +
    facet_wrap(~ in_data) +
    geom_point()
}

ecoregions <- filter(ecoregions, in_data | n_neigh > 0)

if(file.exists('data/ecoregions/poly-nbs-global.rds')) {
  nbs <- readRDS('data/ecoregions/poly-nbs-global.rds')
} else {
  nbs <-
    spdep::poly2nb(pl = ecoregions,
                   row.names = ecoregions$poly_id, # so that names match IDs
                   queen = TRUE) # 1 point in common is sufficient
  length(nbs) == nrow(ecoregions)
  names(nbs) <- attr(nbs, 'region.id') # names do not need to match indices
  
  # add missing neighbors because of polygons that cross the edge of the map
  c('poly 8225', 'poly 8038', 'poly 7987', 'poly 8226', 'poly 7988', 'poly 8040') %in%
    ecoregions$poly_id
  
  layout(matrix(1:6, ncol = 2))
  # 7987 and 8038 are already neighbors
  # 7988 and 8040 are already neighbors
  # 7987 and 8226 do not touch
  add_nb(p1 = 'poly 8225', p2 = 'poly 8226', add = TRUE)
  add_nb(p1 = 'poly 8038', p2 = 'poly 8040', add = TRUE)
  add_nb(p1 = 'poly 7987', p2 = 'poly 7988', add = TRUE)
  add_nb(p1 = 'poly 7987', p2 = 'poly 8040', add = TRUE)
  add_nb(p1 = 'poly 7988', p2 = 'poly 8038', add = TRUE)
  add_nb('poly 8980', 'poly 8987', add = TRUE)
  layout(1)
  saveRDS(nbs, 'data/ecoregions/poly-nbs-global.rds')
}

all(attr(nbs, 'region.id') == ecoregions$poly_id)

d <- mutate(d, poly_id = factor(poly_id, levels = names(nbs)))

# fit the model for estimating mean NDVI ----
if(length(list.files('models/global-test', 'hbam-mean-ndvi-*')) > 1) {
  DATE <- '2025-04-XX'
  m <- readRDS(paste0('models/global-test/hbam-mean-ndvi-', DATE, '.rds'))
} else {
  DATE <- format(Sys.time(), '%Y-%m%-%d-%H-%M-%S')
  
  m <- bam(
    ndvi_15_day_mean ~
      wwf_ecoregion +
      s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
      s(doy, by = wwf_ecoregion, bs = 'cc', k = 10) +
      s(year, by = wwf_ecoregion, bs = 'cr', k = 10) +
      ti(doy, year, by = wwf_ecoregion, bs = c('cc', 'cr'), k = c(5, 10)) +
      s(elevation_m, bs = 'cr', k = 5),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    drop.unused.levels = FALSE,
    discrete = TRUE,
    control = gam.control(nthreads = 1, trace = TRUE))
  
  saveRDS(m, paste0('models/global-test/hbam-mean-ndvi-', DATE, '.rds'))
}

p_hbam <- draw(m, rug = FALSE, - which(smooths(m) == 's(poly_id)'))

ggsave(paste0('figures/hbam-mean-ndvi-', DATE, '.png'), plot = p_hbam,
       width = 8, height = 12, units = 'in', dpi = 600, bg = 'white')

# find pixel-level variance ----
# add fitted values and residuals to the dataset
d$mu_hat <- fitted(m)
d$e <- residuals(m)

# check pixel-level mean residuals
d %>%
  slice_sample(prop = 0.1) %>%
  group_by(long, lat) %>%
  summarize(mean_e = mean(e, na.rm = TRUE), .groups = 'drop') %>%
  pull(mean_e) %>%
  hist()

#' **IF REMOVING PIXEL-LEVEL MEAN RESIDUAL**
if(FALSE) {
  d$e <- d$e - d$mean_e
  d$mu_hat <- d$mu_hat + d$mean_e
}

# calculate pixel-level squared residual (mean e^2 is the variance)
d$e_2 <- d$e^2

saveRDS(d, 'data/hbam-var-ndvi-data.rds')
