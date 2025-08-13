# specs of the EME Linux:
# 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
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
source('analysis/figures/000-default-ggplot-theme.R')

# import water features
# not dropping seas and oceans to remove problematic coast cells
linear <- read_sf('data/water-shapefiles/linear-water.shp')
bodies <- read_sf('data/water-shapefiles/water-polygons.shp')

# import ecoregions
ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp')

#' looks like the widest parts of rivers are included in `bodies`, while
#' `linear` only include minor rivers and canals
filter(linear, grepl('Nile', NAME)) %>%
  st_geometry() %>%
  plot(col = 'blue', axes = TRUE)
filter(bodies, grepl('Nile', Name1)) %>%
  st_geometry() %>%
  plot(col = 'red', border = 'red', add = TRUE)
linear %>%
  st_geometry() %>%
  st_intersection(filter(linear, grepl('Nile', NAME)) %>%
                    st_bbox() %>%
                    st_as_sfc() %>%
                    st_as_sf()) %>%
  plot(col = 'cornflowerblue', add = TRUE)
bodies %>%
  filter(Name1 %in% c('Lake Nasser', 'Nahr an Nil') |
           OBJECTID == '19948') %>%
  st_geometry() %>%
  plot(col = 'darkorange', border = 'darkorange', add = TRUE)

# large bodies of water are split into multiple polygons
plot(st_geometry(filter(bodies, Shape_Area > 7.5)), axes = TRUE,
     col = 'cornflowerblue')

#' not sure what units `bodies$Shape_Area` is in: not km^2 or mi^2
# attempting to see metadata for bodies gives "Unable to transform metadata for item 'e750071279bf450cbd510454a80f2e63' (Error 500)"
# attempting to see metadata for linear gives "Unable to transform metadata for item '273980c20bc74f94ac96c7892ec15aff' (Error 500)"
filter(bodies, grepl('Azov', Name1)) %>%
  mutate(.,
         area_sum = sum(Shape_Area),
         area_km2 = sum(st_area(st_make_valid(.)) / 1e3^2))

# some bodies only have 2 vertices, so they give errors when plotting
if(FALSE) plot(st_geometry(filter(bodies, Shape_Leng == 0)))

head(bodies) %>%
  mutate(n_vertices = map_int(1:n(), \(i) nrow(st_coordinates(.[i, ]))))

filter(bodies, Shape_Leng == 0) %>%
  mutate(n_vertices = map_int(1:n(), \(i) nrow(st_coordinates(.[i, ]))))

# no bodies with negative length
filter(bodies, Shape_Area < 0)

# drop problematic polygons
bodies <- filter(bodies, Shape_Leng > 0)

# some polygons still have invalid geometries
filter(bodies, ISO_CC == 'CA') %>%
  slice(59) %>%
  st_coordinates() %>%
  data.frame() %>%
  ggplot() +
  geom_path(aes(X, Y))

filter(bodies, ISO_CC == 'CA') %>%
  slice(59) %>%
  st_make_valid() %>%
  st_area()

# find polygons with invalid geometries ----
# runs in ~ 21 minutes
tictoc::tic()
bodies <- bodies %>%
  mutate(valid_geom = map_lgl(1:n(), \(i) {
    ! is.na(sf:::CPL_geos_is_valid(st_geometry(.[i, ])))
  }, .progress = TRUE))
tictoc::toc()

# not many polygons with invalid geometry
sum(! bodies$valid_geom)

filter(bodies, ! valid_geom) %>%
  st_make_valid() %>%
  st_area() %>%
  as.numeric() %>%
  `/`(1e6) %>% # convert to km^2
  hist(., xlab = expression(Polygon~area~(km^2)), breaks = 50, main ='')

filter(bodies, ! valid_geom) %>%
  st_make_valid() %>%
  st_geometry() %>%
  plot()

# make invalid polygons valid
ig <- which(! bodies$valid_geom)
bodies[ig, ] <- st_make_valid(bodies[ig, ])
rm(ig)
bodies <- bodies %>%
  mutate(fixed_geom = ! valid_geom) %>%
  select(! valid_geom)

# calculate a raster of proportion of cell covered by water ----
# not masking to keep to find coastlines
r_bodies <- bodies %>%
  vect() %>%
  rasterize(list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/',
                       pattern = '.nc', full.names = TRUE)[1] %>%
              rast(lyr = 'QA'),
            cover = TRUE, background = 0)

plot(r_bodies)

writeRaster(r_bodies, 'data/water-body-raster.tif')
plot(rast('data/water-body-raster.tif'))

# check that masking by prop water and it as a term helps the model
d <- readRDS('data/taiga-test/taiga-ndvi-t-2-s-2-aggr.rds') %>%
  mutate(prop_water = extract(r_bodies, select(., x, y))[, 2])

# check distribution of values to find a good cutoff
plot_grid(
  ggplot(d, aes(prop_water)) +
    geom_histogram(color = 'black', fill = 'grey',
                   breaks = seq(0, 1, by = 0.025), center = 0.0125) +
    labs(x = 'Proportion of pixel covered by water', y = 'Count') +
    scale_x_continuous(expand = c(0.02, 0)),
  ggplot(d, aes(prop_water)) +
    stat_ecdf(geom = 'step', n = 41) +
    geom_vline(xintercept = c(0.4, 0.5), color = 'grey', lty = 'dashed') +
    labs(x = 'Proportion of pixel covered by water', y = 'ECDF') +
    scale_x_continuous(expand = c(0.02, 0)),
  labels = 'AUTO', ncol = 1)
ggsave('figures/taiga-test/prop-water-distr.png', width = 6, height = 8,
       dpi = 600, bg = 'white')

# using 0.4 because 0.5 results in excessive bias, and data loss is minimal
mean(d$prop_water > 0.4) # 0.03519987
d <- filter(d, prop_water < 0.4)

# fits in ~ 5 minutes
m <- bam(
  ndvi ~
    s(y, x, bs = 'sos', k = 200) +
    s(elev_m, bs = 'cr', k = 5) +
    s(year, bs = 'cr', k = 10) +
    s(doy, bs = 'cc', k = 10) +
    s(prop_water, bs = 'cr', k = 5),
  family = gaussian(),
  knots = list(doy = c(0.5, 366.5)),
  data = d,
  method = 'fREML',
  discrete = TRUE,
  nthreads = 10,
  control = gam.control(trace = TRUE))

draw(m, rug = FALSE, nrow = 3)
ggsave('figures/taiga-test/taiga-ndvi-gaussian-sos-aggr-terms-prop-water.png',
       width = 9, height = 9, units = 'in', dpi = 300, bg = 'white')

# make maps of mean and var ----
taiga <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(ecoregion == 'Northwest Territories taiga' |
           ecoregion == 'Muskwa-Slave Lake forests') %>%
  st_union() %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_transform(crs(rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc')))
plot(taiga, col = 'forestgreen')

elevs <- d %>%
  filter(central_date == first(central_date)) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
  rast() %>% # convert RasterLayer to SpatRast
  project(crs(taiga)) %>%
  mask(taiga)

preds_mu <-
  elevs %>%
  as.data.frame(xy = TRUE) %>%
  rename(elev_m = 3) %>%
  filter(! is.na(elev_m)) %>%
  mutate(elev_m = if_else(elev_m < 0, 0, elev_m),
         year = 0, doy = 0, prop_water = 0) %>%
  mutate(mu = predict(object = m, newdata = ., type = 'response',
                      se.fit = FALSE, discrete = FALSE, # T gives low res
                      terms = c('(Intercept)', 's(y,x)', 's(elev_m)')))

preds_s2 <- d %>%
  transmute(x, y, e2 = resid(m)^2) %>%
  summarise(s2 = mean(e2), .by = c(x, y)) %>%
  rast() %>%
  `crs<-`('EPSG:4326') %>%
  project(elevs, res = res(elevs)) %>% #' project to same as `elevs`
  as.data.frame(xy = TRUE)

mean(preds_s2$s2 > 0.04) * 100

plot_grid(
  ggplot(preds_mu) +
    geom_raster(aes(x, y, fill = mu)) +
    geom_sf(data = taiga, fill = 'transparent', color = 'black') +
    scale_fill_viridis_c('NDVI', option = 'A') +
    labs(x = NULL, y = NULL),
  preds_s2 %>%
    mutate(s2 = if_else(s2 > 0.04, 0.04, s2)) %>%
    filter(st_as_sf(., coords = c('x', 'y')) %>%
             st_set_crs('EPSG:4326') %>%
             st_intersects(st_transform(taiga, 'EPSG:4326'),
                           sparse = TRUE) %>%
             map_lgl(\(x) length(x) > 0)) %>%
    ggplot() +
    geom_raster(aes(x, y, fill = s2)) +
    geom_sf(data = taiga, fill = 'transparent', color = 'black') +
    scale_fill_viridis_c(expression(bold(s^'2')), limits = c(0, NA)) +
    labs(x = NULL, y = NULL),
  ncol = 2, labels = 'AUTO')

ggsave('figures/taiga-test/prop-water-estimates.png',
       width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')
