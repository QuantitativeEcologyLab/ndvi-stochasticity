setwd('H/GitHub/ndvi-stochasticity/') # for EME Linux
library('dplyr')         # for data wrangling
library('tidyr')         # for data wrangling
library('sf')            # for shapefiles
library('terra')         # for rasters
library('purrr')         # for functional programming
library('furrr')         # for parallelized functional programming
library('lubridate')     # for working with dates
library('elevatr')       # for elevation data
library('qs')            #'faster than `saveRDS()`/`readRDS()`
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/is_flagged.R')

# using a factor of 2 causes the palearctic group to fail because the dataset
# has too many rows, resulting in the error:
#' `Long vectors are not yet supported. Requested output size must be less than 2147483647.`
s_res <- 2 # factor for spatial aggregation
t_res <- 2 # factor for temporal aggregation

FIRST_100_ONLY <- FALSE # set to TRUE for creating a subset of a few days

if(FIRST_100_ONLY) {
  file_names <-
    list.files('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1',
               pattern = '.nc', full.names = TRUE)[1:100]
  length(file_names)
} else {
  file_names <-
    c(list.files(path = 'data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006',
                 pattern = '.nc', full.names = TRUE, recursive = TRUE),
      list.files(path = 'data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001',
                 pattern = '.nc', full.names = TRUE, recursive = TRUE))
  
  # some years have two rasters
  length(file_names) / 365 # approximate number of years of data
}

# all AVHRR rasters have the same CRS, VIIRS use a different CRS ----
# check crs for first, last, and a random sample of rasters
if(FALSE) {
  file_names %>%
    `[`(if(FIRST_100_ONLY) {
      1:100
    } else {
      c(1, length(file_names), sample(length(file_names), size = 100))
    }) %>%
    map_chr(function(.fn) { # fast enough that it is not worth parallelizing
      crs(rast(.fn))
    }) %>%
    unique() %>%
    cat()
}

# antarctica has too little data and variation to do anything interesting
group_shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  filter(group != 'Antarctic')
if(FALSE) plot(group_shp['group'])

# raster of proportion water (aggregated later, after masking NDVI rasters)
prop_water <- rast('data/water-body-raster.tif') %>%
  ifel(. < 0.4, ., NA_real_)
res(prop_water) # res is 0.05 x 0.05
plot(prop_water)

# download raster of elevation (m) ----
if(file.exists('data/elev-raster.tif'))  {
  elev_m <- rast('data/elev-raster.tif')
} else {
  elev_m <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc') %>%
    get_elev_raster(z = 4) %>% # nearest finer res than 0.05 x 0.05
    rast() %>% # convert from raster to SpatRaster
    project(prop_water, res = res(prop_water)) # convert to same res
  plot(elev_m)
  res(elev_m)
  writeRaster(elev_m, 'data/elev-raster.tif')
}

all(res(prop_water) == res(elev_m))

# calculate slope and aspect (degrees) from the elevation raster ----
if(file.exists('data/slope-raster-degrees.tif')) {
  slope <- rast('data/slope-raster-degrees.tif')
} else {
  slope <- terrain(x = elev_m, v = 'slope', neighbors = 8, unit = 'degrees',
                   filename = 'data/slope-raster-degrees.tif')
  plot(slope)
}

if(file.exists('data/aspect-raster-degrees.tif')) {
  aspect <- rast('data/aspect-raster-degrees.tif')
} else {
  aspect <- terrain(x = elev_m, v = 'aspect', neighbors = 8, unit = 'degrees',
                    filename = 'data/aspect-raster-degrees.tif')
  plot(aspect)
}

# calculate long-term average precipitation raster (1979-01 to 2025-07) ----
# https://cds.climate.copernicus.eu/datasets/ecv-for-climate-change
if(file.exists('data/precipitation-yearly-mean-m-day.tif')) {
  yearly_precip <- rast('data/precipitation-yearly-mean-m-day.tif')
} else {
  r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc')
  
  yearly_precip <-
    list.files(path = 'data/ecmwf-era5-precipitation-rasters',
               pattern = '.grib', full.names = TRUE) %>%
    rast() %>%
    mean() %>%
    project(r_0, res = res(r_0))
  
  writeRaster(yearly_precip, 'data/precipitation-yearly-mean-m-day.tif')
  
  monthly_precip <-
    tibble(fn = list.files(path = 'data/ecmwf-era5-precipitation-rasters',
                           pattern = '.grib', full.names = TRUE),
           r = map(fn, rast),
           month = map_int(r, \(.r) lubridate::month(time(.r)))) %>%
    summarize(mean = list(mean(rast(r))), .by = month) %>%
    pull(mean) %>%
    rast() %>%
    project(r_0, res = res(r_0))
  names(monthly_precip) <- month.name
  
  writeRaster(monthly_precip, 'data/precipitation-monthly-mean-m-day.tif')
}

yearly_precip <- yearly_precip %>%
  mask(st_transform(group_shp)) %>%
  project(prop_water)

# find number of cells per complete raster (i.e., assuming no NAs) ----
# find area of pixels at different latitudes
# lat: 1 deg = ~110 km
# long: 1 deg = 111.320*cos(lat) km
tibble(
  lat = seq(-90, 90, by = 30),
  long = 0,
  lat_deg_to_km = 110,
  long_deg_to_km = 111 * cospi(lat / 180),
  pixel_area_km2 = lat_deg_to_km * long_deg_to_km * 0.05^2)

# create a data frame of all dates
if(file.exists('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds') &
   ! FIRST_100_ONLY) {
  dates <- readRDS('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
} else {
  if(FIRST_100_ONLY) {
    DATES <- as.Date('1981-06-24') + 0:99
  } else {
    DATES <- seq(as.Date('1981-06-24'), as.Date('2025-06-29'), by = 1)
  }
  
  plan(strategy = multisession,
       workers = availableCores(logical = FALSE) - 2)
  
  dates <-
    tibble(date = DATES,
           file_name = future_map_chr(date, function(.date) {
             i <- which(
               # AVHRR
               grepl(format(.date, 'AVH13C1\\.A%Y%j\\.'),file_names)|
                 # VIIRS
                 grepl(format(.date, '_%Y%m%d_'), file_names))
             
             .fn <- file_names[i]
             
             if(length(.fn) == 1) {
               return(.fn)
             } else if(length(.fn) == 2) {
               # sometimes there's 2 satellites & 2 rasters
               return(.fn[2]) # the 2nd is from the more recent satellite
             } else if(length(.fn) == 0) {
               return(NA_character_)
             } else {
               .msg <- paste0('Found ', length(.fn), ' files for ',
                              as.character(.date), '!')
               warning(.msg)
               return(.msg)
             }
           }),
           n_cells = future_map_int(file_name, \(.fn) {
             if(file.exists(.fn)) {
               r <- rast(.fn)
               
               # drop probably/confidently cloudy pixels
               # bits 1 and 0 are 10 or 11
               r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI)
               
               not.na(r$NDVI) %>%
                 values() %>%
                 sum() %>%
                 as.integer() %>%
                 return()
             } else return(NA_integer_)
           }, .progress = TRUE))
  plan(sequential)
  
  if(! FIRST_100_ONLY) {
    saveRDS(dates, 'data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
  }
}

# ensure only 1 raster per date (output should be 0)
sum(grepl('Found', dates$file_name))

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# calculating number of cells and dataset size
n_cells <- sum(dates$n_cells, na.rm = TRUE) # total cells across rasters
max_rows <- 2^31 - 1 # max number of rows for a data frame in R

# largest area is ~7.6 times over max data frame size
data.frame(group_shp) %>%
  summarize(area_prop = sum(area_km2), .by = group) %>%
  mutate(area_prop = area_prop / sum(area_prop),
         rel_dataset_size = area_prop * n_cells / max_rows) %>%
  as_tibble()

# ensure temporal thinning covers all days in a 10 year-period:
all(0:364 %in% ((1:3653 - 1:3653 %% t_res) %% 365))

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

# modeling each group separately
df_sizes <-
  group_shp %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarize(area_1e6_km2 = sum(area_km2) / 1e6) %>%
  mutate(prop_area = area_1e6_km2 / sum(area_1e6_km2),
         # data frame size (n rows), relative to max size: should be < max
         nrow_1e6 = prop_area * (n_cells / s_res^2 / t_res) / 1e6,
         prop_max = nrow_1e6 * 1e6 / max_rows,
         below_max = nrow_1e6 < max_rows / 1e6)
df_sizes

# check number of t_res-day windows with less than t_res days
dates <-
  mutate(dates,
         julian = julian(date),
         date_group = as.integer(julian - (julian %% t_res))) %>%
  # dates and number of rasters for each date_group
  group_by(date_group) %>%
  mutate(n_rasters = sum(! is.na(file_name)),
         start_date = min(date),
         central_date = mean(date),
         end_date = max(date),
         doy = yday(central_date) %>% as.integer(),
         year = year(central_date) %>% as.integer()) %>%
  ungroup()

dates %>%
  group_by(date_group) %>%
  summarise(n_dates = n(),
            n_rasters = unique(n_rasters)) %>%
  group_by(n_dates, n_rasters) %>%
  summarize(total = n())

# make a figure showing the data density (drops 12 groups with no rasters)
p_n_rasters <-
  ggplot(dates) +
  coord_equal(ratio = 4) +
  geom_rect(aes(xmin = doy - 2.5, xmax = doy + 2.5,
                ymin = year - 0.5, ymax = year + 0.5,
                fill = factor(n_rasters))) +
  scale_x_continuous('Day of year', expand = c(0, 0)) +
  scale_y_continuous('Year', expand = c(0, 0)) +
  scale_fill_manual('Number of rasters', values = rev(color('bamako')(5)))+
  theme(legend.position = 'top')

if(! file.exists('figures/input-data/n-rasters-time.pdf')) {
  ggsave('figures/input-data/n-rasters-time.pdf', p_n_rasters,
         width = 9, height = 6, units = 'in', bg = 'white')
}

# create the aggregated datasets ----
# spatRast objects cannot be run in parallel and moved across sessions:
# https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818

# testing object size using one of the datasets
if(FALSE) {
  z0 <- readRDS('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/group-level-datasets/aggregated-indo-malay-oceania-australasia-t-2-s-2-ndvi-data-100-days.rds')
  z1 <- data.frame(z0)
  z2 <- mutate(z0,
               ndvi_aggr = as.integer(ndvi_aggr * 1e4),
               doy = as.integer(doy),
               year = as.integer(year),
               prop_water = as.integer(prop_water * 1e4),
               elev_m = as.integer(elev_m * 1e4),
               slope_deg = as.integer(slope_deg * 1e4),
               aspect_deg = as.integer(aspect_deg * 1e4),
               precip_m_day = as.integer(precip_m_day * 1e4))
  z3 <- select(z2, ! c(date_group, slope_deg, aspect_deg))
  
  # using integers and dropping columns cuts dataset size in half
  tibble(object = paste0('z', 0:3),
         size_Gb = map_dbl(object, \(.n) get(.n) %>%
                             object.size() %>%
                             format(units = 'Gb', digits = 3) %>%
                             gsub(' Gb', '', .) %>%
                             as.numeric()),
         rel_size = size_Gb / size_Gb[1])
}

# ensure all rasters have the same and correct resolution
unique(c(0.05, res(prop_water), res(elev_m), res(slope), res(aspect),
         res(yearly_precip)))

# create a stack of rasters to exract from
# converting to int and dropping terrain variables to minimize dataset size
covariates <-
  list(
    prop_water = aggregate(prop_water, s_res, na.rm=TRUE),
    elev_m = aggregate(elev_m, s_res, na.rm = TRUE) %>%
      ifel(. < 0, 0, .) %>%
      as.int(), # to avoid issues with coast
    # slope_deg = as.int(aggregate(slope, s_res, na.rm = TRUE)),
    # aspect_deg = as.int(aggregate(aspect, s_res, na.rm = TRUE)),
    precip_m_day = aggregate(yearly_precip, s_res, na.rm = TRUE)) %>%
  rast() %>%
  mask(st_transform(group_shp, crs(.)))
plot(covariates)

tibble(x = -119.3949 + rnorm(5), y = 49.93928 + rnorm(5)) %>%
  st_as_sf(crs = crs(covariates), coords = c('x', 'y')) %>%
  extract(covariates, ., ID = FALSE)

# covariates as a data frame to use across cores
covariates_df <- as.data.frame(covariates, xy = TRUE)

# calculate aggregated mean NDVI for each realm ----
GROUPS <- arrange(df_sizes, nrow_1e6)$group

# remove unnecessary RAM-demanding objects
rm(aspect, elev_m, p_n_rasters, prop_water, slope, yearly_precip)

gc()

# for more info on cutoffs, see https://github.com/QuantitativeEcologyLab/stochasticity_parks/blob/main/Scripts/005-create-dataset.R
doy_cutoffs <- tibble(
  date = as.Date(c('2025-01-15', '2025-02-15', '2025-03-15', '2025-04-15',
                   '2025-05-15', '2025-06-15', '2025-07-15', '2025-08-15',
                   '2025-09-15', '2025-10-15', '2025-11-15', '2025-12-15')),
  doy = yday(date),
  month = month(date),
  y = c(60, 65, 70, rep(NA, 5), 80, 69, 62, 59))

# create the datasets
plan(multisession, workers = availableCores(logical = FALSE) - 2)

map_chr(GROUPS, function(.group) {
  # get the shapefile for the area
  .shp <- filter(group_shp, group == .group) %>%
    st_geometry() %>%
    st_as_sf() %>%
    st_transform(crs(rast(file_names[1])))
  
  if(FIRST_100_ONLY) {
    # specify that it's only 100 days of data in the file name
    .filename <-
      paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
             gsub(', ', '-', tolower(.group)),
             '-t-', t_res, '-s-', s_res, '-ndvi-data-100-days.rds')
  } else {
    # create the standard name for a complete dataset
    .filename <-
      paste0('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
             gsub(', ', '-', tolower(.group)),
             '-t-', t_res, '-s-', s_res, '-ndvi-data.qs')
  }
  
  dates %>%
    filter(! is.na(file_name)) %>% # drop missing rasters
    nest(cluster = ! c(date_group, central_date)) %>%
    mutate(
      # import, aggregate temporally, and convert to data frame
      ndvi_df = future_map(cluster, function(.cl) {
        #' import `t_res` rasters for the cluster
        #' but clean each raster individually before averaging
        r_ndvi <-
          map2(.cl$file_name, .cl$date, \(.fn, .date) {
            .r <- rast(.fn)
            
            # drop cloudy pixels
            # not dropping pixels with cloud shadow because it reduces the
            # sample size too much with little change in NDVI: see taiga test
            .r$NDVI <- ifel(is_flagged(.r$QA, 1), NA, .r$NDVI)
            
            # drop pixels with prop_water > 0.4
            .r$NDVI <- # import so we can parallellize
              ifel(project(rast('data/water-body-raster.tif'),
                           .r$NDVI, res = res(.r$NDVI)) < 0.4, .r$NDVI, NA)
            
            # remove unrealistically high NDVI values at high latitudes
            cutoff <- doy_cutoffs$y[which(doy_cutoffs$month == month(.date))]
            
            if(! is.na(cutoff)) {
              .r$NDVI <- ifel(.r$NDVI > 0 & init(.r, 'y') >= cutoff, NA, .r$NDVI)
            }
            
            # aggregating spatially
            .r <- terra::aggregate(.r, s_res, na.rm = TRUE)
            
            # crop and mask based on the shapefile of the group
            .r <- crop(.r, st_transform(.shp, crs(.r)), mask = TRUE)
            
            return(.r$NDVI)
          }) %>%
          rast() %>% # convert list to stack of rasters
          mean(na.rm = TRUE) # aggregate the stack temporally
        
        out <- rast(covariates_df, crs = 'EPSG:4326') %>%
          project(r_ndvi)
        out$NDVI <- r_ndvi
        
        as.data.frame(out, xy = TRUE) %>%
          filter(! is.na(NDVI)) %>%
          return()
      }, .options = furrr_options(seed = NULL), .progress = TRUE)) %>%
    select(date_group, central_date, ndvi_df) %>%
    unnest(ndvi_df) %>%
    rename(ndvi_aggr = NDVI) %>%
    mutate(doy = as.integer(yday(central_date)),
           year = as.integer(year(central_date))) %>%
    qsave(file = .filename, nthreads = availableCores(logical = FALSE) - 2)
  
  return(paste(.group, 'saved.'))
})

plan(sequential)

if(FALSE) { # for testing
  readRDS('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-neotropic-nearctic-t-2-s-2-ndvi-data-100-days.rds') %>%
    filter(central_date <= central_date[1] + t_res * 3) %>%
    ggplot() +
    facet_wrap(~ central_date, ncol = 2) +
    geom_sf(data = read_sf('data/ecoregions/groups-polygons.shp') %>%
              filter(group == 'Neotropic, Nearctic')) +
    geom_raster(aes(x, y, fill = ndvi_aggr)) +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))
}
