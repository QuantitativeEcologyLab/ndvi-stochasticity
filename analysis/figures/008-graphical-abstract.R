library('sf')        # for working with shapefiles
library('terra')     # for working with rasters
library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('lubridate') # for working with dates
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy multi-panel plots
source('analysis/figures/000-default-ggplot-theme.R')

# shapefile of land masses 
shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  filter(group != 'Antarctic') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_union() %>%
  st_as_sf()

files <- list.files(paste0('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/',
                           'raster-files/AVHRR-Land_v006/N07_AVH13C1'),
                    full.names = TRUE)[2:5]

dates <- gsub('.*AVH13C1.A', '', files) %>%
  gsub('.006.*', '', .) %>%
  as_date(format = '%Y%j')

d_ndvi <- purrr::map2(files, dates, function(.file, .date) {
  rast(.file, lyr = 'NDVI') %>%
    mask(shp) %>%
    as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    mutate(date = as_date(.date)) %>%
    as_tibble()
}) %>%
  bind_rows() %>%
  rename(ndvi = NDVI)

#' *SUBSTITUTE RASTERS OF LONG-TERM MEAN AND VARIANCE*
d_sum <- summarize(d_ndvi,
                   mu = mean(ndvi, na.rm = TRUE),
                   s2 = var(ndvi, na.rm = TRUE),
                   .by = c(x, y)) %>%
  mutate(s2 = if_else(s2 > 0.05, 0.05, s2))

p_a <-
  ggplot() +
  facet_wrap(~ date, nrow = 2) +
  geom_sf(data = shp, fill = 'grey', color = 'black') +
  geom_raster(aes(x, y, fill = ndvi), d_ndvi) +
  scale_x_continuous(NULL, breaks = c(), expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = c(), expand = c(0, 0)) +
  scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1)) +
  theme(legend.position = 'top', legend.key.width = rel(2))

p_b <-
  ggplot() +
  geom_sf(data = shp, fill = 'grey') +
  geom_raster(aes(x, y, fill = mu), d_sum) +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_gradientn('Mean NDVI', colors = ndvi_pal, limits = c(-1, 1)) +
  theme(legend.position = 'right')

p_c <-
  ggplot() +
  geom_sf(data = shp, fill = 'grey') +
  geom_raster(aes(x, y, fill = s2), d_sum) +
  scale_x_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(NULL, breaks = NULL, expand = c(0, 0)) +
  scale_fill_davos(name = 'DENVar', limits = c(0, NA), reverse = TRUE) +
  theme(legend.position = 'right')

p <- plot_grid(p_a, plot_grid(p_b, p_c, labels = c('B', 'C'), ncol = 1),
               labels = c('A', ''), nrow = 1, rel_widths = c(1, 1.2))

H <- 531 * 5
W <- 1328 * 5

ggsave('figures/graphical-abstract.png', plot = p, width = W, height = H,
       dpi = 600, units = 'px', bg = 'white')
ggsave('figures/graphical-abstract.pdf', plot = p, width = W, height = H,
       units = 'px', bg = 'white')
