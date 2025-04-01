library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('sf')    # for shapefiles
library('terra') # for rasters
library('purrr') # for functional programming

file_names <-
  list.files(path = 'data/avhrr-viirs-ndvi/raster-files/',
             pattern = '.nc', full.names = TRUE, recursive = FALSE)
length(file_names) / 365 # number of years
sum(grepl('prelim', file_names)) # ensure no rasters are preliminary

# check a few rasters
if(FALSE) {
  plot(rast(file_names[[1]]))
  plot(rast(file_names[[2]])) # compare NDVI to QA
  plot(rast(file_names[[3]]))
  plot(rast(file_names[[4]]))
  plot(rast(file_names[[100]]))
  plot(rast(file_names[[200]]))
} # looks like filtering based on QA has already been done

# no clear connection between QA and NDVI; maybe data are already clean
if(FALSE) {
  r <- rast(file_names[2])
  plot(values(r$QA), values(r$NDVI))
  rm(r)
}

#' one AVHRR raster has ~6M cells, which means only about 330 global rasters fit
#' in a data frame, since the max number of rows is (2^31 - 1) = 2e9
ecoregions <- st_read('data/world-ecosystems/data/commondata/data0/tnc_terr_ecoregions.shp') %>%
  filter(WWF_REALM2 != 'Antarctic') # drop Antarctica

# all rasters have the same CRS
# check crs for first, last, and a random sample of rasters
c(file_names[1], sample(file_names, size = 100), file_names[length(file_names)]) %>%
  map_chr(function(.fn) { # fast enough that it is not worth parallelizing
    crs(rast(.fn), proj = TRUE)
  }) %>%
  unique()

ecoregions <- st_transform(ecoregions, crs(rast(file_names[1]))) %>%
  st_geometry() %>%
  vect()

# create a data frame of all dates
dates <-
  tibble(date = seq(as.Date('1981-06-24'), as.Date('2025-03-15'), by = 1),
         julian = julian(date),
         group = julian - (julian %% 15),
         file_name = map_chr(date, function(.date) {
           .fn <- file_names[grepl(format(.date, '_%Y%m%d_'), file_names)]
           
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
         }))

# ensure no more than two groups (i.e., start and end) have 15 days
dates %>%
  group_by(group) %>%
  summarise(n = n()) %>%
  group_by(n) %>%
  summarize(total = n())

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# need to aggregate temporally to reduce file size before saving ----
# spatRast objects cannot be serialized (i.e., run in parallel): https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818
# data for 2020 omly ----
tictoc::tic() # ~8 minutes on EME Linux
dates_aggr_2020 <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  filter(lubridate::year(date) == 2020) %>%
  nest(cluster = ! group) %>%
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
  }, .progress = 'Aggregating rasters and converting to data frame')) %>%
  select(group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_15_day_mean = mean)

saveRDS(dates_aggr_2020, 'data/ndvi-global-15-day-average-2020-only.rds')
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

# data for 2021 omly ----
tictoc::tic() # ~6 minutes on EME Linux
dates_aggr_2021 <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  filter(lubridate::year(date) == 2021) %>%
  nest(cluster = ! group) %>%
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
  }, .progress = 'Aggregating rasters and converting to data frame')) %>%
  select(group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_15_day_mean = mean)

saveRDS(dates_aggr_2021, 'data/ndvi-global-15-day-average-2021-only.rds')
tictoc::toc()


# all data (~7.55 hours on EME Linux) ----
tictoc::tic()
dates_aggr <-
  dates %>%
  filter(! is.na(file_name)) %>% # drop missing rasters
  nest(cluster = ! group) %>%
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
  select(group, central_date, ndvi_rast) %>%
  unnest(ndvi_rast) %>%
  rename(ndvi_15_day_mean = mean)

saveRDS(dates_aggr, 'data/ndvi-global-15-day-average.rds')
tictoc::toc()

# check the first raster ----
if(FALSE) {
  library('ggplot2')
  source('analysis/figures/000-default-ggplot-theme.R')
  
  dates_aggr %>%
    filter(group == group[1]) %>%
    ggplot() +
    coord_equal() +
    geom_raster(aes(x, y, fill = ndvi_15_day_mean)) +
    scale_fill_gradientn('15-day mean NDVI', colors = ndvi_pal,
                         limits = c(-1, 1)) +
    theme(legend.position = 'top')
}
