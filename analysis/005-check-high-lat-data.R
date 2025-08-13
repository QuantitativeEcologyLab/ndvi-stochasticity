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

# there are some bands of extremely high NDVI for high latitudes
NCORES <- min(availableCores() - 4, 60)
plan(multisession, workers = NCORES)

#' *using the full dataset results in the same, repeated error:*
#' Error in `mutate()`:
#' ℹ In argument: `month_data = future_map(...)`.
#' Caused by error:
#' ℹ In index: 1.
#' Caused by error in `vec_rbind()`:
#' ! Negative `n` in `compact_rep()`.
#' ℹ In file 'utils.c' at line 897.
#' ℹ This is an internal error that was detected in the vctrs package.
if(file.exists('data/avhrr-viirs-ndvi/global-diagnostics/monthly-ecdf-summary-without-ecdf-10-percent-subset.rds')) {
  d_ecdf <- readRDS('data/avhrr-viirs-ndvi/global-diagnostics/monthly-ecdf-summary-without-ecdf-10-percent-subset.rds')
} else {
  d_ecdf <-
    tibble(
      fn = c(
        # AVHRR files
        list.files(
          path = 'data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006',
          recursive = TRUE, full.names = TRUE),
        # VIIRS files
        list.files(
          path = 'data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001',
          recursive = TRUE, full.names = TRUE)),
      date =
        if_else(grepl('AVHRR', fn),
                gsub('.*\\N.._AVH13C1.A', '', fn) %>%
                  gsub('\\..*', '', .) %>%
                  as.Date(format = '%Y%j'),
                substr(fn,
                       nchar(fn) - nchar('YYYYMMDD_cYYYYMMDDHHMMSS.nc') + 1,
                       nchar(fn) - nchar('_cYYYYMMDDHHMMSS.nc')) %>%
                  as.Date(format = '%Y%m%d')),
      year = year(date),
      month = factor(month.name[month(date)], levels = month.name)) %>%
    slice_sample(prop = 0.1) %>%
    nest(month_data = ! month) %>%
    mutate(month_data = future_map(month_data, \(.d) {
      map(.d$fn, \(.fn) {
        .r <- rast(.fn) %>%
          # crop and mask to land area only
          crop(st_transform(ecoregions, crs(.)), mask = TRUE)
        
        # drop cloudy pixels
        .r <- ifel(is_flagged(.r$QA, 1), NA, .r$NDVI)
        
        return(as.data.frame(.r, xy = TRUE))
      }) %>%
        bind_rows() %>%
        summarize(cdf = list(ecdf(NDVI)), .by = y) %>%
        mutate(probs = map(cdf, \(.fun) {
          tibble(ndvi = seq(-0.1, 1, by = 0.01),
                 p = .fun(ndvi))
        })) %>%
        return()
    }, .progress = TRUE, .options = furrr_options(seed = NULL))) %>%
    unnest(month_data) %>%
    unnest(probs)
  
  saveRDS(d_ecdf, 'data/avhrr-viirs-ndvi/global-diagnostics/monthly-ecdf-summary-10-percent-subset.rds')
  
  # saving a version without ECDF to be able to push it to GitHub
  d_ecdf %>%
    select(! cdf) %>%
    saveRDS('data/avhrr-viirs-ndvi/global-diagnostics/monthly-ecdf-summary-without-ecdf-10-percent-subset.rds')
}

# find NDVI values closest to the 99th percentile for each latitude
q_99 <- d_ecdf %>%
  slice(which.min(abs(p - 0.99)), .by = c(month, y)) %>%
  # find NDVI cutoffs after dropping decimals in latitude
  mutate(month = factor(month.name[month], levels = month.name),
         y_2 = floor(y)) %>%
  mutate(ndvi_2 = min(ndvi), .by = c(y_2, month)) %>%
  #' arrange by latitude to make sure `geom_path()` looks right
  arrange(y)

ggplot(q_99) +
  facet_wrap(~ month) +
  geom_path(aes(ndvi_2, y)) +
  scale_x_continuous('NDVI', expand = c(0, 0)) +
  theme(panel.grid = element_blank(), legend.position = 'top',
        text = element_text(face = 'bold'))

cutoffs <- expand_grid(month = unique(q_99$month),
                       y = unique(q_99$y),
                       ndvi = c(0, 1)) %>%
  filter((month == 'January' & ! y < 60) |
           (month == 'February' & ! y < 65) |
           (month == 'March' & ! y < 70) |
           (month == 'September' & ! y < 80) |
           (month == 'October' & ! y < 69) |
           (month == 'November' & ! y < 62) |
           (month == 'December' & ! y < 59)) %>%
  summarize(ymin = min(y), ymax = max(y), .by = c(month, ndvi))

# filter months based on q_99 (black lines or april ref line)
# bright bands simply won't have much data, but the bands are thin anyway
p_cutoffs <-
  ggplot() +
  facet_wrap(~ month) +
  geom_raster(aes(ndvi, y, fill = p), d_ecdf) +
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = ymin, ymax = Inf), cutoffs,
            fill = '#8B102080', color = '#8B0000', linetype = 'dashed') +
  geom_path(aes(ndvi, y), q_99, color = 'white', linewidth = 0.1) +
  geom_path(aes(ndvi_2, y), q_99, color = 'white', linewidth = 0.3) +
  scale_fill_tokyo(name = 'ECDF(NDVI)', reverse = TRUE) +
  scale_x_continuous('NDVI', expand = c(0, 0)) +
  scale_y_continuous('Latitude', expand = c(0, 0), limits = c(-90, 90),
                     breaks = seq(-90, 90, by = 30)) +
  theme(panel.grid = element_blank(), legend.position = 'top',
        text = element_text(face = 'bold'),
        panel.background = element_rect(fill = 'grey'),
        plot.margin = unit(c(5.5, 11, 5.5, 5.5), 'points'))

colorblindr::cvd_grid(p_cutoffs, severity = 1)

ggsave('figures/global-diagnostics/monthly-latitudinal-ecdf.png', p_cutoffs,
       width = 12, height = 10, units = 'in', dpi = 300, bg = 'white')

# the change in cutoff throughout the year can be estimated by a function
# with very little uncertainty, but using these smooth changes risks
# leaving some extreme values in the data. it may be worth doing the check
# at intervals finer than a month (e.g., 7 or 14 days), but we can probably
# make the change in future versions, if necessary.
if(FALSE) {
  doy_cutoffs <- tibble(
    date = c('2025-01-15', '2025-02-15', '2025-03-15', '2025-04-15',
             '2025-05-15', '2025-06-15', '2025-07-15', '2025-08-15',
             '2025-09-15', '2025-10-15', '2025-11-15', '2025-12-15') %>%
      as.Date(),
    doy = yday(date),
    month = month(date),
    y = c(60, 65, 70, rep(NA, 5), 80, 69, 62, 59))
  
  m_cutoffs <- gam(y ~ s(doy, bs = 'cc', k = 7), data = doy_cutoffs,
                   method = 'REML', knots = list(doy = c(0.5, 365.5)))
  plot.gam(m_cutoffs, residuals = TRUE, scheme = 1, pch = 19,
           xlim = c(0.5, 365.5), seWithMean = TRUE,
           trans = \(x) x + coef(m_cutoffs)['(Intercept)'])
  rm(m_cutoffs)
}
