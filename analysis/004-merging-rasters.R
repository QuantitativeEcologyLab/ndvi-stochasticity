library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('lubridate') # for working with dates
source('analysis/figures/000-default-ggplot-theme.R')

file_names <-
  list.files(path = 'data/avhrr-viirs-ndvi/raster-files/',
             pattern = '.nc', full.names = TRUE, recursive = FALSE)
length(file_names) / 365 # approximate number of years of data

# check a few rasters
if(FALSE) {
  plot(rast(file_names[[1]]))
  plot(rast(file_names[[2]])) # compare NDVI to QA
  plot(rast(file_names[[3]]))
  plot(rast(file_names[[4]]))
  plot(rast(file_names[[100]]))
  plot(rast(file_names[[200]]))
  
  # looks like filtering based on QA has already been done
  rast(file_names[[1]]) %>%
    as.data.frame() %>%
    slice_sample(n = 1e4) %>%
    plot(is.na(NDVI) ~ QA, .)
  
  # no clear connection between QA and NDVI; maybe data are already clean
  r <- rast(file_names[2])
  plot(values(r$QA), values(r$NDVI))
  rm(r)
}

# all rasters have the same CRS:
# check crs for first, last, and a random sample of rasters
file_names[c(1, length(file_names),
             sample(length(file_names), size = 100))] %>%
  map_chr(function(.fn) { # fast enough that it is not worth parallelizing
    crs(rast(.fn))
  }) %>%
  unique()

#' `ecoregions` uses same projection as first raster, `file_names[1]`
ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# find number of cells per complete raster (i.e., assuming no NAs)
max_size <- 2^31 - 1 # max data frame size in R
n_cells <-
  rast(file_names[1], lyr = 'QA') %>% # cells in a raster w no NAs
  mask(ecoregions) %>%
  not.na() %>%
  values() %>%
  sum()
n_cells / 1e6 # million land cells per raster (including NAs) 

# lat: 1 deg = ~110 km
# long: 1 deg = 111.320*cos(lat) km
tibble(
  lat = seq(-90, 90, by = 30),
  long = 0,
  lat_deg_to_km = 110,
  long_deg_to_km = 111 * cospi(lat / 180),
  pixel_area_km2 = lat_deg_to_km * long_deg_to_km * 0.05)

# one AVHRR raster has a max of ~8.7M cells
(2^31 - 1) / n_cells # max number of complete rasters per data frame

# find max resolution for datasets to fit in max data frame size
s_res <- 3 # spatial resolution
t_res <- 4 # temporal resolution
((n_cells / t_res^2) * (length(file_names) / s_res)) / max_size

# create a data frame of all dates
if(file.exists('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')) {
  dates <- readRDS('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
} else {
  plan(strategy = multisession,
       workers = availableCores(logical = FALSE) - 2)
  plan()
  dates <-
    tibble(date = seq(as.Date('1981-06-24'), as.Date('2025-05-07'), by =1),
           file_name = future_map_chr(date, function(.date) {
             .fn <- file_names[grepl(format(.date, '_%Y%m%d_'),file_names)]
             
             if(length(.fn) == 1) {
               return(.fn)
             } else if(length(.fn) == 0) {
               return(NA_character_)
             } else {
               .msg <- paste0('Found ', length(.fn), ' files for ',
                              as.character(.date), '!')
               warning(.msg)
               return(.msg)
             }
             return()
           }, .progress = TRUE),
           n_cells = future_map_int(file_name, \(.fn) {
             if(file.exists(.fn)) {
               rast(.fn) %>%
                 mask(ecoregions) %>%
                 not.na() %>%
                 values() %>%
                 sum() %>%
                 as.integer() %>%
                 return()
             } else return(NA_integer_)
           }, .progress = TRUE))
  plan(sequential)
  saveRDS(dates, 'data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
}

sum(grepl('Found', dates$file_name)) # ensure only 1 raster per date

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# calculating number of cells and dataset size
n_cells <- sum(dates$n_cells, na.rm = TRUE) # total cells across rasters
max_rows <- 2^31 - 1 # max number of rows for a data frame in R
n_cells / max_rows # ~164 times over max data frame size
s_res <- 4 # spatial resolution
t_res <- 4 # temporal resolution
n_cells / s_res^2 / t_res / max_rows # modeling as a single dataset

# non-integer time interval gives alternating repetition:
table((0:365 - 0:365 %% t_res))

# check cell sizes before and after spatial aggregation
tibble(
  lat = seq(-90, 90, by = 30),
  long = 0,
  original_lat_km = 110 * 0.05, # original res is 0.05 * 0.05 degrees
  original_long_km = 111 * cospi(lat / 180) * 0.05,
  original_area_km2 = original_lat_km * original_long_km,
  aggr_lat_km = original_lat_km * s_res,
  aggr_long_km = original_long_km * s_res,
  aggr_area_km2 = original_area_km2 * s_res^2)

# modeling each WWF realm separately
df_sizes <-
  ecoregions %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarize(area_1e6_km2 = sum(area_km2) / 1e6) %>%
  mutate(prop_area = area_1e6_km2 / sum(area_1e6_km2),
         # data frame size (n rows), relative to max size: should be < max
         nrow_1e6 = prop_area * (n_cells / s_res^2 / t_res) / 1e6,
         prop_max = nrow_1e6 * 1e6 / max_rows,
         below_max = nrow_1e6 < max_rows / 1e6)
df_sizes

ggplot() +
  geom_sf(data = ecoregions, aes(fill = group)) +
  scale_fill_bright()
  
#' ensure at most two groups (start and end) have < `floor(t_res)` days
dates <-
  mutate(dates,
         julian = julian(date),
         date_group = julian - (julian %% t_res)) %>%
  # dates and number of rasters for each date_group
  group_by(date_group) %>%
  mutate(n_rasters = sum(! is.na(file_name)),
         start_date = min(date),
         central_date = mean(date),
         end_date = max(date),
         doy = yday(central_date),
         year = year(central_date)) %>%
  ungroup()

dates %>%
  group_by(date_group) %>%
  summarise(n_dates = n(),
            n_rasters = unique(n_rasters)) %>%
  group_by(n_dates, n_rasters) %>%
  summarize(total = n())

# make a figure showing the data density (drops 12 groups with no rasters)
dates %>%
  ggplot() +
  coord_equal(ratio = 4) +
  geom_rect(aes(xmin = doy - 2.5, xmax = doy + 2.5,
                ymin = year - 0.5, ymax = year + 0.5,
                # fill = n_rasters)) +
                fill = factor(n_rasters))) +
  scale_x_continuous('Day of year', expand = c(0, 0)) +
  scale_y_continuous('Year', expand = c(0, 0)) +
  scale_fill_manual(
    'Number of rasters',
    values = c('#FAD3BE', '#7A9F9E', '#3572A3', '#21327F', '#190C65')) +
  theme(legend.position = 'top')
ggsave('figures/input-data/n-rasters-time.png',
       width = 8.5, height = 5, units = 'in', dpi = 300, bg = 'white')

# need to aggregate temporally to reduce file size before saving ----
# spatRast objects cannot be serialized (i.e., run in parallel):
# https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818

# calculate aggregated mean NDVI for each realm ---
# cannot serialize rasters across cores, but can serialize realm names
realms <- arrange(df_sizes, nrow_1e6)$WWF_REALM

plan(multisession(workers = min(availableCores() - 2, length(realms))))

future_map_chr(realms, function(.realm) {
  shp <- filter(ecoregions, WWF_REALM == .realm) %>%
    st_geometry() %>%
    st_as_sf()
  
  dates %>%
    filter(! is.na(file_name)) %>% # drop missing rasters
    nest(cluster = ! date_group) %>%
    mutate(
      # aggregate temporally
      ndvi_rast = map(cluster, function(.cl) {
        map(.cl$file_name, \(.fn) rast(.fn, lyr = 'NDVI')) %>% # import
          rast() %>% # convert list to stack of rasters
          mask(shp) %>% # only keep the realm of interest
          mean(na.rm = TRUE) %>% # aggregate temporally
          return()
      }, .progress = 'Calculating mean raster across time'),
      # aggregating spatially
      ndvi_rast = map(ndvi_rast, \(x) {
        as.data.frame(terra::aggregate(x, s_res, na.rm = TRUE), xy = TRUE)
      }, .progress = 'Aggregating spatially and converting to data frame')) %>%
    select(date_group, central_date, ndvi_rast) %>%
    unnest(ndvi_rast) %>%
    rename(ndvi_aggr = mean) %>%
    saveRDS(paste0('data/avhrr-viirs-ndvi/aggregated-', .realm,
                   '-t-', t_res, '-s-', s_res, '-ndvi-data.rds'))
  
  return(paste('Realm', .realm, 'saved.'))
})

plan(sequential)

# data for 2020 only ----
tictoc::tic() # ~8 minutes on EME Linux
dates_aggr_2020 <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  filter(lubridate::year(date) == 2020) %>%
  nest(cluster = ! date_group) %>%
  mutate(central_date = map_dbl(cluster, function(.cl) {
    mean(.cl$date)
  }, .progress = 'Getting mean dates') %>%
    as.Date(), # convert numeric back to date
  ndvi_rast = map(cluster, function(.cl) {
    map(.cl$file_name, \(.fn) {
      rast(.fn, lyr = 'NDVI')
    }) %>%
      rast() %>%
      mask(ecoregions) %>% # drop antarctica and unnecessary pixels
      mean(na.rm = TRUE) %>%
      return()
  }, .progress = 'Calculating mean raster'),
  n_rasters = map_int(cluster, nrow), # find number of rasters per group
  ndvi_rast = map(ndvi_rast, \(x) {
    # aggregating to ensure the final number of rows is below the max
    terra::aggregate(x, s_res) %>%
      as.data.frame(xy = TRUE)
  }, .progress = 'Aggregating rasters and converting to data frame')) %>%
  select(date_group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_aggr_mean = mean)

saveRDS(dates_aggr_2020, 'data/ndvi-global-aggregated-2020-only.rds')
tictoc::toc()

if(FALSE) {
  # check the spatial distribution of the points
  range(dates_aggr_2020$x)
  range(dates_aggr_2020$y)
  
  # good coverage worldwide
  dates_aggr_2020 %>%
    slice_sample(n = 1e5) %>%
    plot(y ~ x, .)
  
  # many days don't have much data on the north pole
  dates_aggr_2020 %>%
    filter(central_date == first(central_date)) %>%
    slice_sample(n = 1e5) %>%
    plot(y ~ x, .)
}

# data for 2021 only ----
tictoc::tic() # ~6 minutes on EME Linux
dates_aggr_2021 <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  filter(lubridate::year(date) == 2021) %>%
  nest(cluster = ! date_group) %>%
  mutate(central_date = map_dbl(cluster, function(.cl) {
    mean(.cl$date)
  }, .progress = 'Getting mean dates') %>%
    as.Date(), # convert numeric back to date
  ndvi_rast = map(cluster, function(.cl) {
    map(.cl$file_name, \(.fn) {
      rast(.fn, lyr = 'NDVI')
    }) %>%
      rast() %>%
      mask(ecoregions) %>% # drop antarctica and unnecessary pixels
      mean(na.rm = TRUE) %>%
      return()
  }, .progress = 'Calculating mean raster'),
  n_rasters = map_int(cluster, nrow), # find number of rasters per group
  ndvi_rast = map(ndvi_rast, \(x) {
    # aggregating to ensure the final number of rows is below the max
    terra::aggregate(x, s_res) %>%
      as.data.frame(xy = TRUE)
  }, .progress = 'Aggregating rasters and converting to data frame')) %>%
  select(date_group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_aggr_mean = mean)

saveRDS(dates_aggr_2021, 'data/ndvi-global-aggregated-2021-only.rds')
tictoc::toc()


# all data (~7.55 hours on EME Linux using 15-day aggregation) ----
tictoc::tic()
dates_aggr <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  nest(cluster = ! date_group) %>%
  mutate(central_date = map_dbl(cluster, function(.cl) {
    mean(.cl$date)
  }, .progress = 'Getting mean dates') %>%
    as.Date(), # convert numeric back to date
  ndvi_rast = map(cluster, function(.cl) {
    map(.cl$file_name, \(.fn) {
      rast(.fn, lyr = 'NDVI')
    }) %>%
      rast() %>%
      mask(ecoregions) %>% # drop antarctica and unnecessary pixels
      mean(na.rm = TRUE) %>%
      return()
  }, .progress = 'Calculating mean raster'),
  n_rasters = map_int(cluster, nrow), # find number of rasters per group
  ndvi_rast = map(ndvi_rast, \(x) {
    # aggregating to ensure the final number of rows is below the max
    terra::aggregate(x, 2) %>%
      as.data.frame(xy = TRUE)
  }, .progress = 'Aggregating rasters')) %>%
  select(date_group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_aggr_mean = mean)

saveRDS(dates_aggr, 'data/ndvi-global-aggregated.rds')
tictoc::toc()

# check the first raster ----
if(FALSE) {
  library('ggplot2')
  source('analysis/figures/000-default-ggplot-theme.R')
  
  dates_aggr %>%
    filter(date_group == date_group[1]) %>%
    ggplot() +
    coord_equal() +
    geom_raster(aes(x, y, fill = ndvi_aggr_mean)) +
    scale_fill_gradientn('Aggregated mean NDVI', colors = ndvi_pal,
                         limits = c(-1, 1)) +
    theme(legend.position = 'top')
}
