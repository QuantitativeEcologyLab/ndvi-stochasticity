library('terra')  # for rasters
library('sf')     # for simple features
library('ggplot2') # for fancy plots
library('dplyr')  # for data wrangling
source('analysis/figures/000-default-ggplot-theme.R')
r <- rast("H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810701_c20170609165932.nc", lyr = 'NDVI')
eco <- st_read('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(WWF_MHTNAM != 'Inland Water') %>%
  filter(ECO_NAME != 'Great Lakes') %>%
  st_geometry()

plot(r)
plot(eco, add = TRUE)

d <- as.data.frame(r, xy = TRUE)
nrow(d) / 1e6

p <- ggplot() +
  geom_sf(data = eco, color = 'black', fill = 'grey30', lwd = 0.1) +
  geom_raster(aes(x, y, fill = NDVI), d) +
  geom_sf(data = eco, color = 'black', fill = 'transparent', lwd = 0.1) +
  scale_x_continuous(NULL, expand = c(0, 0)) +
  ylab(NULL) +
  scale_fill_gradientn(colours = ndvi_pal, limits = c(-1, 1)) +
  theme(legend.position = 'top')

ggsave('figures/ndvi-1981-07-01-raw-raster.png',
       p, width = 8, height = 4, units = 'in', bg = 'white', dpi = 300)
