# create list of cell neighbors
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
ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(group != 'antarctica') %>%
  st_geometry() %>%
  st_transform(crs(r_0)) %>%
  st_as_sf()

#' *WILL LIKELY NEED TO SPLIT BY GROUPS*

# all rasters use same coords (about 1.6 M cells with values; 5.3 M total)
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc',
            lyr = 'QA') %>% # to have a value for all cell regardless of raster
  aggregate(2) # because aggregated in the dataset
values(r_0) <- as.integer(cells(r_0)) # add cell id

r_0 <- mask(crop(r_0, ecoregions), ecoregions)
names(r_0) <- 'cell_id'
plot(r_0)
writeRaster(r_0, 'data/cell-id.tif')

plot(rast('data/cell-id.tif'))

# create a list of neighbors
tictoc::tic() # run time on EME Linux: XXX hours
nbs <- nbs_from_rast(r_0)
head(nbs)
all(range(values(r_0), na.rm = TRUE) == range(as.numeric(names(nbs))))
saveRDS(nbs, 'data/cell-nbs-list.rds')
tictoc::toc()

# run a quick test to check names
test_sf <- slice(ecoregions, 100) %>%
  st_geometry()
test_r <- crop(rast('data/global-cell-nbs.tif'), test_sf)
plot(test_r)
plot(test_sf, add = TRUE, lwd = 2, col = '#00000040')
min(values(test_r), na.rm = TRUE)
nbs[as.character(min(values(test), na.rm = TRUE))]
rm(test)
