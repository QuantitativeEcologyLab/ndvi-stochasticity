#' *THE MODELS FIT IN THIS SCRIPT ARE FOR INITIAL TESTING ONLY*
#' *DO NOT USE THESE MODELS TO PREDICT*
library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('sf')      # for shapefiles
library('mgcv')    # for Generalized Additive Models
library('sf')      # for simple feature objects
library('terra')   # for rasters
library('elevatr') # for extracting elevation
library('gratia')  # for plotting GAMs
source('analysis/figures/000-default-ggplot-theme.R')

# shapefile of ecoregions for masking rasters
ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

if(file.exists('data/hbam-ndvi-data-2020-2021-only.rds')) {
  d <- readRDS('data/hbam-ndvi-data-2020-2021-only.rds')
} else {
  # resolution is (0.05 degrees)^2 = (3 degree minutes)^2
  r_0 <- rast('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc',
              lyr = 'QA') %>% # to have a value for each cell irrespective of raster
    aggregate(2) # because aggregated in the dataset
  values(r_0) <- 1
  r_0 <- mask(r_0, st_transform(ecoregions, crs(r_0)))
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
    r_elev <- project(r_elev, r_0) # project to common projection and resolution
    plot(r_elev)
    writeRaster(r_elev, 'data/elev-raster.tif')
  }
  
  d <- bind_rows(readRDS('data/ndvi-global-15-day-average-2020-only.rds'),
                 readRDS('data/ndvi-global-15-day-average-2021-only.rds')) %>%
    #' use `exact_extract()` for more data?
    mutate(doy = lubridate::yday(central_date),
           elevation_m = extract(r_elev, data.frame(x, y))[, 2]) %>%
    na.omit() # drop some NAs (< 0.1%)
  saveRDS(d, 'data/hbam-ndvi-data-2020-2021-only.rds')
}

# set excessively low elevations to sea level (occurs mostly near islands)
range(d$elevation_m)
if(FALSE) {
  plot(st_geometry(ecoregions))
  d %>%
    filter(elevation_m < 0) %>%
    points(y ~ x, ., col = 'red', pch = '.', cex = 2)
}

d <- mutate(d, elevation_m = if_else(elevation_m < 0, 0, elevation_m))

nbs <- readRDS('data/cell-nbs-list.rds')

d <-
  mutate(d,
         cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs)))

# global smooths only
# 2020 only: initial is 34 s, fit is 5 s
# 2020 and 2021: initial is 40 s, fit is 10 s
# 2020 and 2021 on EME linux: initial is 50 s, fit is 10 s
if(file.exists('models/global-test/m0-2020-2021-gam.rds')) {
  m0 <- readRDS('models/global-test/m0-2020-2021-gam.rds')
} else {
  m0 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(doy, bs = 'cc', k = 10),
            family = gaussian(), # much faster than betar() with similar results
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            nthreads = 10,
            control = gam.control(trace = TRUE))
  draw(m0, rug = FALSE)
  saveRDS(m0, 'models/global-test/m0-2020-2021-gam.rds')
}

# adding s(doy, ecoregion); no common doy effect across all ecoregions
# 2020 only: initial is 23 s, fit is 38 s (almost all in first iteration)
# 2020 and 2021: initial is 50 s, fit is 80 s (almost all in first iteration)
# 2020 and 2021 on EME linux using 10 threads: initial is 46 s, fit is 51 s
if(file.exists('models/global-test/m1-2020-2021-gam.rds')) {
  m1 <- readRDS('models/global-test/m1-2020-2021-gam.rds')
} else {
  m1 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10),
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            samfrac = 0.01,
            nthreads = 10,
            control = gam.control(trace = TRUE))
  draw(m1, rug = FALSE)
  saveRDS(m1, 'models/global-test/m1-2020-2021-gam.rds')
}

# adding polygon_id MRF (much, much slower...)
# 2020 only: initial is 53 s, fit is 1.3 hours
# 2020 and 2021: initial is XX s, fit is 1.6 hours
# 2020 and 2021 on EME linux using 10 threads: initial is ~60 s, fit is 20 min
if(file.exists('models/global-test/m2-2020-2021-gam.rds')) {
  m2 <- readRDS('models/global-test/m2-2020-2021-gam.rds')
} else {
  m2 <- bam(ndvi_15_day_mean ~
              s(elevation_m, bs = 'cr', k = 5) +
              s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
              s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            drop.unused.levels = FALSE,
            discrete = TRUE,
            nthreads = 10,
            control = gam.control(trace = TRUE))
  draw(m2, rug = FALSE, select = 1:2)
  saveRDS(m2, 'models/global-test/m2-2020-2021-gam.rds')
}

# unexplained deviance, starting from second iteration (it's 0 at iteration 1)
# for 2020 data only
plot(2:14,
     c(630903, 630326, 629312, 628184, 627524, 627288, 627211, 627188, 627174,
       627171, 627170, 627170, 627170), type = 'l',
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
if(file.exists('models/global-test/m3-2020-2021-gam.rds')) {
  m3 <- readRDS('models/global-test/m3-2020-2021-gam.rds')
} else {
  m3 <- bam(ndvi_15_day_mean ~
              wwf_ecoregion + # to avoid shrinkage of intercepts
              s(elevation_m, bs = 'cr', k = 5) +
              s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
              poly_id, # to check if not using MRF saves a lot of time
            family = gaussian(),
            data = d,
            method = 'fREML',
            knots = list(doy = c(0.5, 366.5)),
            discrete = TRUE,
            nthreads = 10,
            control = gam.control(trace = TRUE))
  
  draw(m3, rug = FALSE)
  saveRDS(m3, 'models/global-test/m3-2020-2021-gam.rds')
}

# compare deviance explained
tibble(model = 0:3,
       dev_expl = map_dbl(list(m0, m1, m2, m3), \(.m) {
         100 - .m$deviance / .m$null.deviance * 100
       }),
       smooths = map_chr(list(m0, m1, m2, m3), function(.m) {
         paste(smooths(.m), collapse = ', ')
       }))

# model dev_expl terms
#   m0    14.3 % s(elevation_m), s(doy)
#   m1    57.5 % s(elevation_m), s(doy,wwf_ecoregion)
#   m2    67.9 % s(elevation_m), s(doy,wwf_ecoregion), s(poly_id)
#   m3    67.9 % s(elevation_m), s(doy):wwf_ecoregion + wwf_ecoregion + poly_id

# test for variance -----
# add fitted values and residuals to the dataset
d$mu_hat <- fitted(m2)
d$e <- residuals(m2)

# check pixel-level mean residuals
d <- d %>%
  group_by(x, y) %>%
  mutate(mean_e = mean(e, na.rm = TRUE)) %>%
  ungroup()

d %>%
  group_by(x, y) %>%
  slice(1) %>%
  pull(mean_e) %>%
  hist(breaks = 50, main = 'Mean pixel-level residual')

#' *IF REMOVING PIXEL-LEVEL MEAN RESIDUAL*
d <- mutate(d,
            e = e - mean_e,
            e_2 = e^2, # squared residuals (mean e^2 is var)
            mu_hat = mu_hat + mean_e)

# fit the model for the variance
m_var <- bam(e_2 ~
               s(elevation_m, bs = 'cr', k = 5) +
               s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
               s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
             family = gaussian(),
             data = d,
             method = 'fREML',
             knots = list(doy = c(0.5, 366.5)),
             drop.unused.levels = FALSE,
             discrete = TRUE,
             nthreads = 10,
             control = gam.control(trace = TRUE))
draw(m_var, rug = FALSE, select = 1:2)
saveRDS(m_var, 'models/global-test/m-var-2020-2021-gam.rds')

# compare MRF to sos smoothers ----
d <- d %>%
  filter(central_date == central_date[1]) %>%
  select(x, y, ndvi_15_day_mean, poly_id)

# k = 1e3: fits in ~ 30 seconds
# k = 2e3: fits in ~ 653 seconds
m_sos <- bam(ndvi_15_day_mean ~ s(y, x, bs = 'sos', k = 2e3), # s(lat,long)
             family = gaussian(),
             data = d,
             method = 'fREML',
             discrete = TRUE,
             control = gam.control(trace = TRUE))
saveRDS(m_sos, 'models/global-test/test-mean-gam-sos-only.rds')

png('figures/test-hgams/test-hgam-sos.png', width = 7.5, height = 10,
    units = 'in', bg = 'white', res = 300)
layout(matrix(c(rep(1, 4), 2:5), ncol = 2, byrow = TRUE))
plot(m_sos, scheme = 4, too.far = 0.01, n2 = 200)
plot(m_sos, theta = 0, phi = 90, n2 = 100, too.far = 0.01)
plot(m_sos, theta = 90, phi = 90, n2 = 100, too.far = 0.01)
plot(m_sos, theta = 180, phi = 90, n2 = 100, too.far = 0.01)
plot(m_sos, theta = 270, phi = 90, n2 = 100, too.far = 0.01)
dev.off()

# fits in XXX seconds
tictoc::tic()
m_mrf <- bam(ndvi_15_day_mean ~ s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
             family = gaussian(),
             data = d,
             method = 'fREML',
             discrete = TRUE,
             control = gam.control(trace = TRUE))
tictoc::toc()

# fits in > tens of minutes (not fit)
tictoc::tic()
m_fe <- bam(ndvi_15_day_mean ~ poly_id,
            family = gaussian(),
            data = d,
            method = 'fREML',
            control = gam.control(trace = TRUE))
tictoc::toc()

# fits in 326 seconds (5.43 minutes), basis size for REs is 3739 / 3951
tictoc::tic()
m_re <- bam(ndvi_15_day_mean ~ s(poly_id, bs = 're'),
            family = gaussian(),
            data = d,
            method = 'fREML',
            discrete = TRUE,
            control = gam.control(trace = TRUE))
tictoc::toc()
(1 - m_re$deviance / m_re$null.deviance) * 100
m_re$df.null - m_re$df.residual

# test for model of the variance
d <- mutate(d, e2 = resid(m_sos)^2)

m_sos_var <- bam(
  e2 ~ s(y, x, bs = 'sos', k = 2e3), # s(lat, long)
  family = gaussian(),
  data = d,
  method = 'fREML',
  discrete = TRUE,
  control = gam.control(trace = TRUE))
saveRDS(m_sos_var, 'models/global-test/test-var-gam-sos-only.rds')
