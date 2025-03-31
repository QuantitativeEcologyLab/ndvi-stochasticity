library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('sf')      # for simple fleature objects
library('terra')   # for rasters
library('elevatr') # for extracting elevation
library('gratia')  # for plotting GAMs
source('analysis/figures/default-ggplot-theme.R')
source('functions/add_nb.R')

ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# resolution is (0.05 degrees)^2 = (3 degree minutes)^2
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc',
            lyr = 'QA') %>% # to have a value for each cell irrespective of raster
  aggregate(2) # because aggregated in the dataset
values(r_0) <- 1
plot(r_0)
res(r_0)

# create raster of ecoregion
if(file.exists('data/ecoregions/wwf-ecoregions.tif')) {
  r_eco <- rast('data/ecoregions/wwf-ecoregions.tif')
} else {
  r_eco <-
    ecoregions %>%
    select(WWF_MHTNAM) %>%
    vect() %>%
    rasterize(r_0, field = 'WWF_MHTNAM')
  plot(r_eco, legend = FALSE)
  writeRaster(r_eco, 'data/ecoregions/wwf-ecoregions.tif')
}

# create raster of polygon ID
if(file.exists('data/ecoregions/ecoregion-polygon-id.tif')) {
  r_poly_id <- rast('data/ecoregions/ecoregion-polygon-id.tif')
} else {
  r_poly_id <-
    ecoregions %>%
    select(poly_id) %>%
    vect() %>%
    rasterize(r_0, field = 'poly_id')
  plot(r_poly_id, legend = FALSE)
  writeRaster(r_poly_id, 'data/ecoregions/ecoregion-polygon-id.tif')
}

# get raster of elevations
if(file.exists('data/elev-raster.tif')) {
  r_elev <- rast('data/elev-raster.tif')
} else {
  r_elev <- get_elev_raster(ecoregions, z = 3) #'`z=3` gives `0.076 < res(r_0)`
  r_elev <- rast(r_elev) #' `{elevatr}` v.0.99.0 still uses `raster::raster()`
  res(r_elev)
  res(r_0)
  plot(r_elev)
  r_elev <- project(r_elev, r_0) # project to common projection and resolution
  writeRaster(r_elev, 'data/elev-raster.tif')
}

# create raster of distance from coast
if(file.exists('data/distance-from-coast-m.tif')) {
  r_dist <- rast('data/distance-from-coast-m.tif')
} else {
  r_dist <- mask(r_0, vect(ecoregions), inverse = TRUE) %>%
    distance()
  r_dist <- mask(r_dist, st_transform(ecoregions, crs(r_dist)))
  plot(r_dist, main = 'Distance from coast (m)')
  plot(log10(r_dist), main = expression(bold(log['10'](Distance~(m)))))
  writeRaster(r_dist, 'data/distance-from-coast-m.tif')
}

if(file.exists('data/hbam-ndvi-data-2020-2021-only.rds')) {
  d <- readRDS('data/hbam-ndvi-data-2020-2021-only.rds')
} else {
  d <- bind_rows(readRDS('data/ndvi-global-15-day-average-2020-only.rds'),
                 readRDS('data/ndvi-global-15-day-average-2021-only.rds')) %>%
    #' use `exact_extract()` for more data?
    mutate(doy = lubridate::yday(central_date),
           poly_id = extract(r_poly_id, data.frame(x, y))[, 2], # 1 is pixel ID
           wwf_ecoregion = extract(r_eco, data.frame(x, y))[, 2],
           distance_coast_m = extract(r_dist, data.frame(x, y))[, 2],
           elevation_m = extract(r_elev, data.frame(x, y))[, 2]) %>%
    # convert strings to factors
    mutate(wwf_ecoregion = factor(wwf_ecoregion)) %>%
    na.omit() # drop some NAs (< 0.1%)
  saveRDS(d, 'data/hbam-ndvi-data-2020-2021-only.rds')
}

# set excessively low elevations to sea level (occurs mostly near islands)
if(FALSE) {
  plot(st_geometry(ecoregions))
  d %>%
    filter(elevation_m < -20) %>%
    points(y ~ x, ., col = 'red', pch = 19)
}

d <- mutate(d, elevation_m = if_else(elevation_m < -20, 0, elevation_m))

# create list of neighbors ----
ecoregions$in_data <- ecoregions$poly_id %in% unique(d$poly_id)
filter(ecoregions, ! in_data) # some polygons have no data
filter(ecoregions, ! in_data & n_neigh == 0) # some polygons have no data and have no neighbors

if(FALSE) {
  ggplot(ecoregions, aes(log10(area_km2), in_data)) + geom_point()
  
  ggplot(ecoregions, aes(log10(area_km2), n_neigh)) +
    facet_wrap(~ in_data) +
    geom_point()
}

ecoregions <- filter(ecoregions, in_data | n_neigh > 0)

nbs <- spdep::poly2nb(pl = ecoregions,
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

all(attr(nbs, 'region.id') == ecoregions$poly_id)

d <- mutate(d, poly_id = factor(poly_id, levels = names(nbs)))

# global smooths only
# 2020 only: initial is 34 s, fit is 5 s
# 2020 and 2021: initial is 40 s, fit is 10 s
if(file.exists('models/global-test/m0-2020-2021.rds')) {
  m0 <- readRDS('models/global-test/m0-2020-2021.rds')
} else {
  m0 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(log(distance_coast_m), bs ='cr', k = 5) +
              s(doy, bs = 'cc', k = 10),
            family = gaussian(), # much faster than betar() with similar results
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            control = gam.control(nthreads = 1, trace = TRUE))
  draw(m0, rug = FALSE)
  saveRDS(m0, 'models/global-test/m0-2020-2021.rds')
}

# adding s(doy, ecoregion); no common doy effect across all ecoregions
# 2020 only: initial is 23 s, fit is 38 s (almost all in first iteration)
# 2020 and 2021: initial is 50 s, fit is 80 s (almost all in first iteration)
if(file.exists('models/global-test/m1-2020-2021-gam.rds')) {
  m1 <- readRDS(m1, 'models/global-test/m1-2020-2021-gam.rds')
} else {
  m1 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(log(distance_coast_m), bs ='cr', k = 5) +
              s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10),
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            control = gam.control(nthreads = 1, trace = TRUE))
  draw(m1, rug = FALSE)
  saveRDS(m1, 'models/global-test/m1-2020-2021-gam.rds')
}

# adding polygon_id MRF (much, much slower...)
# 2020 only: initial is 53 s, fit is 1.3 hours
# 2020 and 2021: initial is 11 s, fit is 1.6 hours
if(file.exists('models/global-test/m2-2020-2021-gam.rds')) {
  m2 <- readRDS(m2, 'models/global-test/m2-2020-2021-gam.rds')
} else {
  m2 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(log(distance_coast_m), bs ='cr', k = 5) +
              s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
              s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            drop.unused.levels = FALSE,
            discrete = TRUE,
            control = gam.control(nthreads = 1, trace = TRUE))
  draw(m2, rug = FALSE, select = -4)
  saveRDS(m2, 'models/global-test/m2-2020-2021-gam.rds')
}

# unexplained deviance, starting from second iteration (it's 0 at iteration 1)
# for 2020 data only
plot(2:16,
     c(350145, 349113, 346434, 342350, 338839, 337245, 336736, 336517,
       336447, 336425, 336418, 336416, 336415, 336415, 336415), type = 'l',
     ylab = 'Deviance', xlab = 'Iteration')

d %>%
  filter(central_date == central_date[1]) %>%
  mutate(.,
         mu = predict(m2, ., type = 'response', se.fit = FALSE)) %>%
  ggplot(aes(x, y, fill = mu)) +
  coord_equal() +
  geom_raster() +
  scale_fill_gradientn('Mean NDVI', colours = ndvi_pal, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = paste('Estimated NDVI on', d$central_date[1]))

# testing ecoregion and polygon as fixed effects to hopefully decrease fitting time
#' 2020 AND 2021: setup is 70 s, fit is 18 minutes
if(file.exists('models/global-test/m3-2020-2021.rds')) {
  m3 <- readRDS('models/global-test/m3-2020-2021.rds')
} else {
  m3 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(log(distance_coast_m), bs ='ad', k = 10) +
              s(doy, by = wwf_ecoregion, bs = 'cc', k = 10) +
              poly_id,
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            control = gam.control(nthreads = 1, trace = TRUE))
  draw(m3, rug = FALSE)
  saveRDS(m3, 'models/global-test/m3-2020-2021.rds')
}

# compare deviance explained
tibble(model = 0:3,
       dev_expl = map_dbl(list(m0, m1, m2, m3), \(.m) {
         100 - .m$deviance / .m$null.deviance * 100
       }))
