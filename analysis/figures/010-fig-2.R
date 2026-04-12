library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('mgcv')    # for predicting from models
library('gratia')  # for working with GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots with multiple panels
source('analysis/figures/000-default-ggplot-theme.R')
source('analysis/figures/000-robinson-objects.R')

r_eco <- read_sf('data/ecoregions/ecoregions-polygons-edited.shp') %>%
  filter(group == 'Neotropic, Nearctic') %>%
  vect() %>%
  rasterize(rast('data/water-body-raster.tif'), field = 'biome')

shp_0 <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  st_geometry() %>%
  st_union() %>%
  st_as_sf()

shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  filter(group != 'Antarctic') %>%
  st_geometry() %>%
  st_union() %>%
  st_as_sf()
plot(shp, col = 'darkgreen')

d_rob <-
  rast('output/long-term-preds.tif') %>%
    mask(st_transform(shp, crs(.))) %>%
    as.data.frame(xy = TRUE) %>%
    as_tibble() %>%
    mutate(prop_water = extract(rast('data/water-body-raster.tif'),
                                tibble(x, y))[, 2],
           elev_m = extract(rast('data/elev-raster.tif'),
                            tibble(x, y))[, 2]) %>%
    filter(prop_water < 0.4)

d_rob <- d %>%
  rast(crs = crs(rast('output/long-term-preds.tif'))) %>%
  project(robinson_crs) %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble() %>%
  mutate(biome = extract(project(r_eco, robinson_crs), tibble(x, y))[, 2])

# Robinson projections of long-term mean NDVI and DENVar
fig_2 <- plot_grid(
  ggplot() +
    geom_sf(data = bounds, fill = 'white', color = 'black') +
    geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
    geom_raster(aes(x, y, fill = mu_hat), d_rob) +
    geom_sf(data = shp_0, color = 'black', fill = 'transparent',
            lwd = 0.05) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn('Mean NDVI', colours = ndvi_pal,
                         limits = c(-1, 1)) +
    theme_void() +
    theme(legend.position = 'top', panel.grid = element_blank(),
          text = element_text(face = 'bold'), legend.key.width = rel(2)),
  ggplot() +
    geom_sf(data = bounds, fill = 'white', color = 'black') +
    geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
    geom_raster(aes(x, y, fill = s2_hat), d_rob) +
    geom_sf(data = shp_0, color = 'black', fill = 'transparent',
            lwd = 0.05) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_davos(name = 'DENVar', limits = c(0, NA), reverse = TRUE) +
    theme_void() +
    theme(legend.position = 'top', panel.grid = element_blank(),
          text = element_text(face = 'bold'), legend.key.width = rel(2)),
  labels = 'AUTO', ncol = 1)

ggsave('figures/fig-2.png', fig_2, width = 8.5, height = 10, dpi = 600,
       bg = 'white')

ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
  geom_raster(aes(x, y, fill = biome == 'Rock and Ice' & ! is.na(biome)),
              d_rob) +
  geom_sf(data = shp_0, color = 'black', fill = 'transparent',
          lwd = 0.05) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual('Rock and Ice', values = c('grey30', 'cornflowerblue'),
                    labels = c('No', 'Yes')) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))
