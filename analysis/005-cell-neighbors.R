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

# all rasters use same coords (about 1.6 M cells with values; 5.3 M total)
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc',
            lyr = 'QA') %>% # to have a value for all cell regardless of raster
  aggregate(2) # because aggregated in the dataset
values(r_0) <- as.integer(cells(r_0))

# shapefile of ecoregions for masking rasters
ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp') %>%
  st_transform(crs(r_0))

r_0 <- mask(crop(r_0, ecoregions), ecoregions)
names(r_0) <- 'cell_id'
plot(r_0)
writeRaster(r_0, 'data/cell-id.tif')

if(FALSE) { # for testing
  ecoregions <- ecoregions[100, ]
  
  # plotting as characters for easier comparison
  # not saving as characters because it replaces values with factor level number
  plot(`values<-`(r_0, as.character(values(r_0))))
}

tictoc::tic() # run time on EME Linux: XXX hours
nbs_df <-
  adjacent(r_0, cells = cells(r_0), directions = 8, include = TRUE) %>%
  as.data.frame() %>%
  transmute(cell_id = values(r_0, lyr = 'cell_id', na.rm = TRUE),
            ref_cell = V1, # first column is the starting cell
            # add the 8 surrounding neighbors
            # not using furrr because cannot parallellize raster operations
            adjacent = map(1:n(), \(i) {
              .z <- c(V2[i], V3[i], V4[i], V5[i], V6[i], V7[i], V8[i], V9[i])
              
              .values <- map_int(.z, \(.cell_id) {
                if(is.nan(.cell_id)) {
                  return(NA)
                } else {
                  return(r_0[.cell_id]$cell_id[1])
                }})

              .values <- .values[! is.na(.values)] # drop cells w NA values
              
              if(length(.values) == 0) {
                return(0)
              } else {
                return(as.character(.values))
              }
            })) %>%
  as_tibble() # for easier viewing when printing
nbs_df
nbs <- nbs_df$adjacent
names(nbs) <- as.character(nbs_df$cell_id)
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

