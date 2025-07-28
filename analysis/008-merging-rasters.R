library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('lubridate') # for working with dates
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/is_flagged.R')

file_names <-
  list.files(path = 'data/avhrr-viirs-ndvi/raster-files/',
             pattern = '.nc', full.names = TRUE, recursive = FALSE)
length(file_names) / 365 # approximate number of years of data

# all rasters have the same CRS ----
# check crs for first, last, and a random sample of rasters
file_names[c(1, length(file_names),
             sample(length(file_names), size = 100))] %>%
  map_chr(function(.fn) { # fast enough that it is not worth parallelizing
    crs(rast(.fn))
  }) %>%
  unique()

#' `ecoregions` uses same projection as first raster, `file_names[1]`
ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp')

# find number of cells per complete raster (i.e., assuming no NAs) ----
# lat: 1 deg = ~110 km
# long: 1 deg = 111.320*cos(lat) km
tibble(
  lat = seq(-90, 90, by = 30),
  long = 0,
  lat_deg_to_km = 110,
  long_deg_to_km = 111 * cospi(lat / 180),
  pixel_area_km2 = lat_deg_to_km * long_deg_to_km * 0.05^2)

# create a data frame of all dates
if(file.exists('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')) {
  dates <- readRDS('data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
} else {
  plan(strategy = multisession,
       workers = availableCores(logical = FALSE) - 2)
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
               r <- rast(.fn)
               
               # drop probably/confidently cloudy pixels
               # bits 1 and 0 are 10 or 11
               r$NDVI <- ifel(is_flagged(r$QA, 1), NA, r$NDVI)
               
               #' *DROP PIXELS WITH PROP_WATER > 0.4*
               
               not.na(r$NDVI) %>%
                 values() %>%
                 sum() %>%
                 as.integer() %>%
                 return()
             } else return(NA_integer_)
           }, .progress = TRUE))
  plan(sequential)
  saveRDS(dates, 'data/avhrr-viirs-ndvi/ndvi-raster-metadata.rds')
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

# largest area is ~5 times over max data frame size
data.frame(ecoregions) %>%
  summarize(area_prop = sum(area_km2), .by = group) %>%
  mutate(area_prop = area_prop / sum(area_prop),
         rel_dataset_size = area_prop * n_cells / max_rows) %>%
  as_tibble()

s_res <- 2 # spatial resolution
t_res <- 2 # temporal resolution

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

# modeling each WWF realm separately
# paleoartctic should be ok because it tends to have more NA values
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

# check number of 2-days windows with less than 2 days
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
  coord_equal(ratio = 3) +
  geom_rect(aes(xmin = doy - 2.5, xmax = doy + 2.5,
                ymin = year - 0.5, ymax = year + 0.5,
                fill = factor(n_rasters))) +
  scale_x_continuous('Day of year', expand = c(0, 0)) +
  scale_y_continuous('Year') +
  scale_fill_manual(
    'Number of rasters',
    values = c('#FFE599', '#5F7D13', '#003F4C')) +
  theme(legend.position = 'top')

if(! file.exists('figures/input-data/n-rasters-time.pdf')) {
  ggsave('figures/input-data/n-rasters-time.pdf',
         width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')
}

# create the aggregated datasets ----
# spatRast objects cannot be run in parallel and moved across sessions:
# https://stackoverflow.com/questions/67445883/terra-package-returns-error-when-try-to-run-parallel-operations/67449818#67449818

# dropping points that are at the edge of boundaries to remove excessive
# variance near coasts. this means we cannot expand the spatial extent
# of each group slightly to reduce boundary issues with GAMs (since we
# can't both expand and contract the extents).
if(FALSE) {
  z <- rast(file_names[3], lyr = 'NDVI')
  shp <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
    filter(realm == 'Nearctic') %>%
    slice(18) %>%
    st_cast('POLYGON', warn = FALSE) %>%
    st_geometry() %>%
    st_as_sf() %>%
    slice(100)
  plot(shp)
  z <- crop(z, shp)
  values(z) <- rnorm(ncell(z))
  plot(z)
  
  layout(matrix(1:4, ncol = 2, byrow = TRUE))
  # plot all cells that are in or touch the polygon
  plot(crop(z, shp, mask = TRUE, touches = TRUE), main = 'touches = TRUE')
  plot(shp, add = TRUE, border = 'red', lwd = 2)
  
  # plot all cells that are in polygon
  plot(crop(z, shp, mask = TRUE, touches = FALSE), main = 'touches = FALSE')
  plot(shp, add = TRUE, border = 'red', lwd = 2)
  
  # plot all cells that are 90% inside the shapefile
  ifel(rasterize(shp, z, cover = TRUE, background = 0) > 0.9, z, NA) %>%
    plot(main = 'Cells 90% inside')
  plot(shp, add = TRUE, border = 'red', lwd = 2)
  
  # plot all cells that are 99.9% inside the shapefile
  ifel(rasterize(shp, z, cover = TRUE, background = 0) >= 0.999, z, NA) %>%
    plot(main = 'Cells 99.9% inside')
  plot(shp, add = TRUE, border = 'red', lwd = 2)
}

# calculate aggregated mean NDVI for each realm ----
GROUPS <- arrange(df_sizes, nrow_1e6)$group

# parallelizing across realm names requires too much RAM
map_chr(GROUPS, function(.group) {
  shp <- filter(ecoregions, group == .group) %>%
    st_geometry() %>%
    st_as_sf()
  
  dates %>%
    filter(! is.na(file_name)) %>% # drop missing rasters
    nest(cluster = ! c(date_group, central_date)) %>%
    mutate(
      # import and aggregate temporally
      ndvi_rast = map(cluster, function(.cl) {
        #' import `t_res` rasters for the cluster
        map2(.cl$file_name, .cl$date, \(.fn, .date) {
          .r <- rast(.fn, lyr = 'NDVI') %>%
            mask(shp, touches = FALSE)
          .month <- month(.date)
          
          # remove unrealistically high NDVI values at high latitudes
          if(.month %in% c(1:4, 9:12)) {
            .lats <- # create a raster with TRUE if above max latitude 
              init(.r, 'y') >= case_when(.month %in% c(1:4, 11, 12) ~ 60,
                                         .month %in% 9:10 ~ 70)
            .r <- ifel(.r > 0.2 & .lats, NA, .r)
          }
          
          # aggregating spatially
          .r <- terra::aggregate(.r, s_res, na.rm = TRUE)
          
          
          return(.r)
        }) %>%
          rast() %>% # convert list to stack of rasters
          mean(na.rm = TRUE) %>% # aggregate temporally
          return()
      }, .progress = 'Calculating mean raster across time and space'),
      ndvi_rast = map(ndvi_rast, \(x) {
        as.data.frame(x, xy = TRUE)
      }, .progress = 'Converting to data frame')) %>%
    select(date_group, central_date, ndvi_rast) %>%
    unnest(ndvi_rast) %>%
    rename(ndvi_aggr = mean) %>%
    saveRDS(paste0('data/avhrr-viirs-ndvi/aggregated-', .group,
                   '-t-', t_res, '-s-', s_res, '-ndvi-data.rds'))
  
  return(paste('Group', .group, 'saved.'))
})

plan(sequential)

if(FALSE) { # for testing
  readRDS('data/avhrr-viirs-ndvi/aggregated-islands-t-2-s-2-ndvi-data.rds') %>%
    filter(central_date <= central_date[1] + t_res * 2) %>%
    ggplot(aes(x, y, fill = ndvi_aggr)) +
    coord_sf(crs = 'EPSG:4326') +
    facet_wrap(~ central_date, ncol = 1) +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))
}

