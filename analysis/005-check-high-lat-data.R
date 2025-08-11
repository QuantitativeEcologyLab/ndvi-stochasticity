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
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R') # get_legend() from cowplot v. 1.1.3 fails
source('functions/nbs_from_rast.R') # gives a list of neighboring cells

ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp')
taiga <- filter(ecoregions, ecoregion == 'Northwest Territories taiga')

# some very large NDVI values at the edge of the preliminary cleaning ----
plot(st_geometry(ecoregions))
plot(st_geometry(taiga), add = TRUE, col = 'red')

list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/',
           'N07_AVH13C1.A198201',
           full.names = TRUE) %>%
  map(\(x) rast(x, lyr = 'NDVI')) %>%
  rast() %>%
  crop(taiga, mask = TRUE) %>%
  plot()

# import raster for 1982-01-01 (AVHRR data)
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A19820101.006.2022271093046.nc') %>%
  crop(st_transform(ecoregions, crs(.)), mask = TRUE)
plot(r_0) # no clear banding

decoded_0 <- decode_qa(as.numeric(values(r_0$QA)), sensor = 'AVHRR')

# make rasters for each QA parameter
for(cn in colnames(decoded_0)[3:ncol(decoded_0)]) {
  r_0[[cn]] <- r_0$QA
  values(r_0[[cn]]) <- decoded_0[[cn]]
  r_0[[cn]] <- mask(r_0[[cn]], r_0$QA)
}
plot(r_0)

# import raster for 2024-01-10 (VIIRS data)
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001/VIIRS-Land_v001_JP113C1_NOAA-20_20240110_c20240124141257.nc') %>%
  crop(st_transform(ecoregions, crs(.)), mask = TRUE)
plot(r_0) # clear banding visible

decoded_0 <- decode_qa(as.numeric(values(r_0$QA)), sensor = 'VIIRS')

# make rasters for each QA parameter
for(cn in colnames(decoded_0)[3:ncol(decoded_0)]) {
  r_0[[cn]] <- r_0$QA
  values(r_0[[cn]]) <- decoded_0[[cn]]
  r_0[[cn]] <- mask(r_0[[cn]], r_0$QA)
}
plot(r_0)

d_0 <- tibble(as.data.frame(r_0, na.rm = TRUE, xy = TRUE))
range(d_0$y)
summary(filter(d_0, NDVI > 0.2)) # max lat is too high for such high NDVI
plot(NDVI ~ y, filter(d_0, y > 60))

# snow/ice flag is not very reliable alone: many values in snow are > 0.25
d_0 %>%
  ggplot() +
  facet_grid(. ~ land_type) +
  geom_boxplot(aes(snow_ice, NDVI), alpha = 0.01)

d_0 %>%
  summarise(n = n(),
            median = median(NDVI, na.rm = TRUE),
            q_0.75 = quantile(NDVI, 0.75, na.rm = TRUE),
            upr_IQR = quantile(NDVI, 0.75) + 1.5 * IQR(NDVI, na.rm = TRUE),
            .by = c(snow_ice))
rm(r_0, d_0)

# first year of data, every 10 days
list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/',
           pattern = '.nc',
           full.names = TRUE)[seq(1, 365, by = 10)] %>%
  map(\(fn) rast(fn, lyr = 'NDVI')) %>%
  rast() %>%
  mask(ecoregions) %>%
  plot()

# make some rasters of summary stats every 10 rasters for each month ----
file_names_10 <- tibble(
  fn = c(
    # AVHRR files
    list.files(
      path = 'data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006',
      # julian day of year ends in 0: AYYYYJJ0 (take one every 10 days)
      pattern = 'A......0', recursive = TRUE, full.names = TRUE),
    # VIIRS files
    list.files(
      path = 'data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001',
      # day of month ends in 1: _YYYYMMD1 (approximately one every 10 days)
      pattern = '_.......1_', recursive = TRUE, full.names = TRUE)),
  date = if_else(grepl('AVHRR', fn),
                 gsub('.*\\N.._AVH13C1.A', '', fn) %>%
                   gsub('\\..*', '', .) %>%
                   as.Date(format = '%Y%j'),
                 substr(fn,
                        nchar(fn) - nchar('YYYYMMDD_cYYYYMMDDHHMMSS.nc') + 1,
                        nchar(fn) - nchar('_cYYYYMMDDHHMMSS.nc')) %>%
                   as.Date(format = '%Y%m%d')),
  year = year(date),
  month = month(date))
range(file_names_10$date) # ensure dates are ok and not NA

# max(NDVI); runs in ~ 30 minutes on EME Linux ----
if(! dir.exists('data/avhrr-viirs-ndvi/global-diagnostics'))
  dir.create('data/avhrr-viirs-ndvi/global-diagnostics')

if(! file.exists('data/avhrr-viirs-ndvi/global-diagnostics/max-raster-jan.tif')) {
  map(1:12, \(i) {
    r_max <-
      file_names_10 %>%
      filter(month == i) %>%
      pull(fn) %>%
      map(function(.fn) mask(rast(.fn, lyr = 'NDVI'), ecoregions)) %>%
      rast() %>%
      max(na.rm = TRUE)
    
    writeRaster(r_max, paste0('data/avhrr-viirs-ndvi/global-diagnostics/',
                              'max-raster-', tolower(month.abb[i]), '.tif'),
                overwrite = TRUE)
  }, .progress = TRUE)
}

# plot the max values
max_values <- tibble(
  month_int = 1:12,
  month_name = factor(month.name[month_int], levels = month.name),
  nested_data = map(month_int, function(i) {
    paste0('data/avhrr-viirs-ndvi/global-diagnostics/max-raster-',
           tolower(month.abb[i]), '.tif') %>%
      rast() %>%
      mask(ecoregions) %>%
      as.data.frame(xy = TRUE, na.rm = TRUE)
  })) %>%
  unnest(nested_data)

# rasters for each month
p_max_rast <-
  ggplot() +
  coord_sf(crs = 'EPSG:4326') +
  facet_wrap(~ month_name) +
  geom_raster(aes(x, y, fill = max), max_values) +
  labs(x = NULL, y = NULL) +
  scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))

ggsave('figures/global-diagnostics/monthly-max-values.png', p_max_rast,
       width = 15, height = 5.5, units = 'in', dpi = 300, bg = 'white')

# histograms of max values at high latitudes
p_max_hist <-
  max_values %>%
  filter(max > 0.2, y > 50) %>%
  ggplot(aes(y)) +
  facet_wrap(~ month_name) +
  geom_histogram(bins = 500) +
  scale_x_continuous('Latitude', limits = c(50, 85)) +
  scale_y_continuous('Cells with max(NDVI) > 0')

ggsave('figures/global-diagnostics/monthly-max-values-above-0.2-hist.png',
       p_max_hist, width = 15, height = 10, units = 'in', dpi = 300,
       bg = 'white')


} else {
}

