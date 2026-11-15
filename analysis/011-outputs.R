# NOTE: predicting at the 0.05-degree resolution despite the models using
# data at the 0.2-degree resolution, since the scale of the input
# data is not the scale of the GAMs' bases
library('sf')      # for shapefiles
library('dplyr')   # for data wrangling
library('mgcv')    # for predicting from models
library('gratia')  # for selecting model terms
library('terra')   # for working with rasters
library('ggplot2') # for fancy plots
library('qs')      # for importing models quickly
library('tidyr')   # for nesting data
library("lubridate")
source('analysis/figures/000-default-ggplot-theme.R')

GROUPS <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  filter(group != 'Antarctic') %>%
  pull(group)

# long-term estimates ----
for (g in GROUPS) {
  cat(paste0('Importing ', g, '...\n'))
  formatted_g <- gsub(',', '', gsub(' ', '-', tolower(g)))
  
  m_mu <- paste0('bam-mean-ndvi-', formatted_g) %>%
    list.files(., path = 'models/global-models', full.names = TRUE) %>%
    qread(nthreads = 40)
  
  m_s2 <- paste0('bam-variance-ndvi-', formatted_g) %>%
    list.files(., path = 'models/global-models', full.names = TRUE) %>%
    qread(nthreads = 40)
  
  shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
    filter(group == g) %>%
    st_geometry() %>%
    st_as_sf()
  
  r_prop_water <- rast('data/water-body-raster.tif') %>%
    ifel(. > 0.4, NA, .) %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE)
  
  r_elev <- rast('data/elev-raster.tif') %>%
    project(r_prop_water) %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE)
  
  plot(r_prop_water)
  plot(r_elev)
  
  EXCLUDE <- smooths(m_mu)[grepl('year', smooths(m_mu)) |
                             grepl('doy', smooths(m_mu))]
  preds <- as.data.frame(r_prop_water, xy = TRUE, na.rm = TRUE) %>%
    as_tibble() %>%
    rename(prop_water = layer) %>%
    mutate(elev_m = extract(r_elev, tibble(x, y))[, 2],
           year = 1e3, doy = 1e3) %>% # irrelevant since excluded
    mutate(mu_hat = predict(m_mu, newdata = ., type = 'response',
                            se.fit = FALSE, discrete = FALSE,
                            exclude = EXCLUDE) %>%
             as.numeric()) %>%
    mutate(s2_hat = predict(m_s2, newdata = ., type = 'response',
                            se.fit = FALSE, discrete = FALSE,
                            exclude = EXCLUDE) %>%
             as.numeric()) %>%
    mutate(year = NA, doy = NA)
  preds
  
  readr::write_csv(preds, paste0('output/long-term-estimates-',
                                 formatted_g,'.csv'))
  
  ggplot(preds, aes(x, y, fill = mu_hat)) +
    coord_sf(crs = 'EPSG:4326') +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
  
  ggplot(preds, aes(x, y, fill = s2_hat)) +
    coord_sf(crs = 'EPSG:4326') +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_davos('DENVar')
}

d <-
  purrr::map(list.files('output', 'long-term-estimates',
                        full.names=TRUE) %>%
               setdiff(., 'output/long-term-estimates.tif'),
             function(.fn) {
               readr::read_csv(.fn, show_col_types = FALSE) %>%
                 select(x, y, mu_hat, s2_hat) %>%
                 rast() %>%
                 terra::`crs<-`(crs(rast('data/elev-raster.tif'))) %>%
                 project(rast('data/elev-raster.tif')) %>%
                 as.data.frame(xy = TRUE)
             }, .progress = TRUE) %>%
  bind_rows() %>%
  filter(y <= 70) %>% # already filtered for prop_water < 0.4
  as_tibble()

d %>%
  rast() %>%
  `crs<-`(crs(rast('data/elev-raster.tif'))) %>%
  writeRaster('output/long-term-estimates.tif')

hist(rast('output/long-term-estimates.tif'))
plot(rast('output/long-term-estimates.tif'))

ggplot(d, aes(x, y, fill = mu_hat)) +
  coord_sf(crs = 'EPSG:4326') +
  geom_raster() +
  labs(x = NULL, y = NULL) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))

ggplot(d, aes(x, y, fill = s2_hat)) +
  coord_sf(crs = 'EPSG:4326') +
  geom_raster() +
  labs(x = NULL, y = NULL) +
  scale_fill_davos(name = 'DENVar', limits = c(NA_real_, 0.05))

# yearly rasters ----
for (g in GROUPS) {
  cat(paste0('Importing ', g, '...\n'))
  formatted_g <- gsub(',', '', gsub(' ', '-', tolower(g)))
  
  m_mu <- paste0('bam-mean-ndvi-', formatted_g) %>%
    list.files(., path = 'models/global-models', full.names = TRUE) %>%
    qread(nthreads = 40)
  
  m_s2 <- paste0('bam-variance-ndvi-', formatted_g) %>%
    list.files(., path = 'models/global-models', full.names = TRUE) %>%
    qread(nthreads = 40)
  gc()
  
  shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
    filter(group == g) %>%
    st_geometry() %>%
    st_as_sf()
  
  r_prop_water <- rast('data/water-body-raster.tif') %>%
    ifel(. > 0.4, NA, .) %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE)
  
  r_elev <- rast('data/elev-raster.tif') %>%
    project(r_prop_water) %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE)
  
  EXCLUDE <- smooths(m_mu)[grepl('doy', smooths(m_mu))]
  preds <- as.data.frame(r_prop_water, xy = TRUE, na.rm = TRUE) %>%
    as_tibble() %>%
    rename(prop_water = layer) %>%
    expand_grid(year = 1981:2025) %>%
    mutate(elev_m = extract(r_elev, tibble(x, y))[, 2],
           doy = 1e3) %>% # irrelevant since excluded
    mutate(mu_hat = predict(m_mu, newdata = ., type = 'response',
                            se.fit = FALSE, discrete = FALSE,
                            exclude = EXCLUDE) %>%
             as.numeric()) %>%
    mutate(s2_hat = predict(m_s2, newdata = ., type = 'response',
                            se.fit = FALSE, discrete = FALSE,
                            exclude = EXCLUDE) %>%
             as.numeric()) %>%
    mutate(doy = NA)
  
  readr::write_csv(preds, paste0('output/yearly-estimates-',
                                 formatted_g,'.csv'))
}

d <-
  purrr::map(list.files('output', 'yearly-estimates',
                        full.names=TRUE) %>%
               setdiff(., 'output/yearly-estimates.csv'),
             function(.fn) {
               readr::read_csv(.fn, show_col_types = FALSE) %>%
                 select(x, y, mu_hat, s2_hat, year) %>%
                 nest(y_data = ! year) %>%
                 mutate(y_data = purrr::map(y_data, \(.d) {
                   rast(.d) %>%
                     terra::`crs<-`(crs(rast('data/elev-raster.tif'))) %>%
                     project(rast('data/elev-raster.tif')) %>%
                     as.data.frame(xy = TRUE) %>%
                     filter(y <= 70) # already filtered for prop_water < 0.4
                 })) %>%
                 unnest(y_data)
             }, .progress = TRUE) %>%
  bind_rows() %>%
  as_tibble()

readr::write_csv(d, 'output/yearly-estimates.csv')
range(readr::read_csv('output/yearly-estimates.csv', num_threads = 10)$year) ==
  c(1981, 2025)

d %>%
  group_by(year) %>%
  nest(y_data = ! year) %>%
  mutate(y_data = purrr::map(y_data, function(.d) {
    .d %>%
      select(x, y, mu_hat) %>%
      rast()
  })) %>%
  pull(y_data) %>%
  `names<-`(1981:2025) %>%
  rast() %>%
  writeRaster("output/yearly-estimates-mean-ndvi.tif", overwrite = TRUE)

d %>%
  group_by(year) %>%
  nest(y_data = ! year) %>%
  mutate(y_data = purrr::map(y_data, function(.d) {
    .d %>%
      select(x, y, s2_hat) %>%
      rast()
  })) %>%
  pull(y_data) %>%
  `names<-`(1981:2025) %>%
  rast() %>%
  writeRaster("output/yearly-estimates-denvar.tif", overwrite = TRUE)

# daily estimates ----
gc()

purrr::map_int(1982:2025, \(x) dir.create(paste0("H:/ndvi-denvar-daily-estimates/", x)))

DATES <- seq(as.Date("1981-01-01"), as.Date("2025-12-31"), by = 1)

# drop dates manually if prediction gets interrupted.
# the code for automatic drop is not worth the time because it requires
# possibly listing too many files.

shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  select(group) %>%
  filter(group != "Antarctic") %>%
  st_as_sf()

r_prop_water <- rast('data/water-body-raster.tif') %>%
  ifel(. > 0.4, NA, .) %>%
  crop(st_transform(shp, crs = crs(.)), mask = TRUE)

r_elev <- rast('data/elev-raster.tif') %>%
  project(r_prop_water) %>%
  crop(st_transform(shp, crs = crs(.)), mask = TRUE)

r_group <- rasterize(shp, r_prop_water, field = "group")
plot(r_group)

newd <- as.data.frame(c(r_prop_water, r_group, r_elev), xy = TRUE, na.rm = TRUE) %>%
  as_tibble() %>%
  rename(prop_water = layer,
         elev_m = file37905ccd4e9d) %>%
  filter(y < 70)

m_mu_a <- qread("models/global-models/bam-mean-ndvi-africa-2026-02-03.qs", nthreads = 40)
m_mu_i <- qread("models/global-models/bam-mean-ndvi-indo-malay-oceania-australasia-2026-05-20.qs", nthreads = 40)
m_mu_n <- qread("models/global-models/bam-mean-ndvi-neotropic-nearctic-2025-10-29.qs", nthreads = 40)
m_mu_p <- qread("models/global-models/bam-mean-ndvi-palearctic-2026-03-24.qs", nthreads = 40)
m_s2_a <- qread("models/global-models/bam-variance-ndvi-africa-2026-02-03.qs", nthreads = 40)
m_s2_i <- qread("models/global-models/bam-variance-ndvi-indo-malay-oceania-australasia-2026-05-20.qs", nthreads = 40)
m_s2_n <- qread("models/global-models/bam-variance-ndvi-neotropic-nearctic-2026-02-24.qs", nthreads = 40)
m_s2_p <- qread("models/global-models/bam-variance-ndvi-palearctic-2026-03-24.qs", nthreads = 40)

map (DATES, function(.date) {
  file_name <- paste0(year(.date), "-", yday(.date), ".tif")
  cat(paste0('Predicting for ', .date, ' (file ', file_name, ')...\n'))
  
  preds <- newd %>%
    mutate(doy = yday(.date),
           year = year(.date)) %>%
    nest(dat = !group) %>%
    mutate(dat = purrr::map2(group, dat, function(.g, .d) {
      m_mu <- get(paste0("m_mu_", tolower(substr(.g, 1, 1))))
      m_s2 <- get(paste0("m_s2_", tolower(substr(.g, 1, 1))))
      dat$mu_hat <- predict.bam(m_mu, newdata = dat, type = "response",
                                se.fit = FALSE, discrete = FALSE) %>%
        as.numeric() # to drop matrix format
      dat$s2_hat <- predict.bam(m_s2, newdata = dat, type = "response",
                                se.fit = FALSE, discrete = FALSE) %>%
        as.numeric()
      dat
    })) %>%
    tidyr::unnest(dat)
  
  preds %>%
    select(x, y, mu_hat, s2_hat) %>%
    rast() %>%
    terra::writeRaster(preds,
                       paste0('H:/ndvi-denvar-daily-estimates/',
                              year(.date), "/",
                              year(.date), "-", yday(.date),'.tif'))
})
