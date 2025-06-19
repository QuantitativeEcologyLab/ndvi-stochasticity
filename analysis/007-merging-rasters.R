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

# all rasters have the same CRS ----
# check crs for first, last, and a random sample of rasters
file_names[c(1, length(file_names),
             sample(length(file_names), size = 100))] %>%
  map_chr(function(.fn) { # fast enough that it is not worth parallelizing
    crs(rast(.fn))
  }) %>%
  unique()

#' `ecoregions` uses same projection as first raster, `file_names[1]`
ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# find number of cells per complete raster (i.e., assuming no NAs) ----
# lat: 1 deg = ~110 km
# long: 1 deg = 111.320*cos(lat) km
tibble(
  lat = seq(-90, 90, by = 30),
  long = 0,
  lat_deg_to_km = 110,
  long_deg_to_km = 111 * cospi(lat / 180),
  pixel_area_km2 = lat_deg_to_km * long_deg_to_km * 0.05)

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

# ensure only 1 raster per date (output should be 0)
sum(grepl('Found', dates$file_name))

# some dates are missing a raster (not available on the server)
all(! is.na(dates$file_name))
mean(! is.na(dates$file_name))
filter(dates, is.na(file_name))

# calculating number of cells and dataset size
n_cells <- sum(dates$n_cells, na.rm = TRUE) # total cells across rasters
max_rows <- 2^31 - 1 # max number of rows for a data frame in R

# lasgest area is ~65 times over max data frame size
data.frame(ecoregions) %>%
  summarize(area_prop = sum(area_km2), .by = group) %>%
  mutate(area_prop = area_prop / sum(area_prop),
         rel_dataset_size = area_prop * n_cells / max_rows) %>%
  as_tibble()

s_res <- 4 # spatial resolution
t_res <- 4 # temporal resolution

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
# preview the groups
if(FALSE) {
  ggplot(ecoregions, aes(fill = group)) +
    geom_sf(lwd = 0.05, col = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_bright(name = 'Group', labels = stringr::str_to_sentence) +
    theme(legend.position = 'top')
  ggsave('figures/input-data/polygon-groups.png',
         width = 10, height = 4.6, units = 'in', dpi = 1200, bg = 'white')
}

GROUPS <- arrange(df_sizes, nrow_1e6)$group

# cannot serialize rasters across cores, but can serialize realm names
plan(multisession(workers = min(availableCores() - 2, length(GROUPS))))

future_map_chr(GROUPS, function(.group) {
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
            mask(shp)
          .month <- month(.date)
          
          # remove unrealistically high NDVI values at high latitudes
          if(.month %in% c(1:4, 9:12)) {
            .lats <- # create a raster with TRUE if above max latitude 
              init(.r, 'y') >= case_when(.month %in% c(1:4, 11, 12) ~ 60,
                                         .month %in% 9:10 ~ 70)
            .r <- ifel(.r > 0.2 & .lats, NA, .r)
          }
          
          return(.r)
        }) %>%
        rast() %>% # convert list to stack of rasters
        # mask(shp) %>% # only keep the realm of interest
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
  saveRDS(paste0('data/avhrr-viirs-ndvi/aggregated-', .group,
                 '-t-', t_res, '-s-', s_res, '-ndvi-data.rds'))

return(paste('Group', .group, 'saved.'))
})

plan(sequential)

if(FALSE) { # for testing
  readRDS('data/avhrr-viirs-ndvi/aggregated-islands-t-4-s-4-ndvi-data.rds') %>%
    filter(central_date <= central_date[1] + t_res * 2) %>%
    ggplot(aes(x, y, fill = ndvi_aggr)) +
    coord_sf(crs = 'EPSG:4326') +
    facet_wrap(~ central_date, ncol = 1) +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = c(-1, 1))
}
