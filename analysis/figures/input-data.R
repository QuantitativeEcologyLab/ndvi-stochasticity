library('dplyr')     # for data wrangling
library('sf')        # for simple features
library('terra')     # for rasters
library('ggplot2')   # for fancy plots
library('tidyterra') # for plotting spatRasters with ggplot2
library('khroma')    # for colorblind-friendly palettes
source('analysis/figures/default-ggplot-theme.R')

eco <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(WWF_MHTNAM != 'Inland Water') %>%
  mutate()

if(FALSE) {
khroma::info() %>%
  filter(type == 'qualitative') %>%
  pull(palette)
  
  # potential discrete palettes
  plot_scheme(color('bright')(7))
  plot_scheme(color('highcontrast')(3))
  plot_scheme(color('vibrant')(7))
  plot_scheme(color('discreterainbow')(16))
  
  khroma::info() %>%
    filter(type == 'sequential') %>%
    pull(palette)
}

# individual polygons ----
n_distinct(eco$poly_id)
pal_polys <- c(color('okabeito')(8), 'grey')
plot_scheme_colorblind(pal_polys)

p_poly <-
  ggplot(eco) +
  geom_sf(fill = pal_polys[1:nrow(eco) %% length(pal_polys) + 1],
          lwd = .05) +
  scale_x_continuous(expand = c(0, 0))

ggsave('figures/input-data/polygons.png', p_poly,
       width = 10, height = 4.4, units = 'in', dpi = 600, bg = 'white')

# ecoregions ----
n_distinct(eco$WWF_MHTNAM)

p_eco <-
  ggplot(eco) +
  geom_sf(aes(fill = WWF_MHTNAM), color = 'black', lwd = .05) +
  scale_fill_discreterainbow(name = 'Ecoregion') +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5))
ggsave('figures/input-data/ecoregions.png', p_eco,
       width = 10.1, height = 6, units = 'in', dpi = 600, bg = 'white')

# elevation ----
r_elev <- rast('data/elev-raster.tif') %>%
  crop(., st_transform(eco, crs(.)), mask = TRUE)

hist(r_elev)
values(r_elev) <- if_else(values(r_elev) < -400, -400, values(r_elev))

ggplot() +
  geom_spatraster(data = r_elev < 0)

p_elev <-
  ggplot() +
  geom_spatraster(data = r_elev) +
  geom_sf(data = eco, color = 'black', fill = 'transparent', lwd = 0.05) +
  scale_fill_lajolla(name = 'Elevation above sea level (m)') +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5))
p_elev
ggsave('figures/input-data/elev-m.png', p_elev,
       width = 10, height = 5, units = 'in', dpi = 600, bg = 'white')

# distance from coast?
r_dist <- rast('data/distance-from-coast-m.tif') %>%
  crop(., st_transform(eco, crs(.)), mask = TRUE)

p_dist <-
  ggplot() +
  geom_spatraster(data = r_dist / 1e3) +
  geom_sf(data = eco, color = 'black', fill = 'transparent', lwd = 0.05) +
  scale_fill_davos(name = 'Distance from coast (km)', reverse = TRUE,
                   limits = c(0, NA), range = c(0.1, 1)) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5))
p_dist

ggsave('figures/input-data/dist-km.png', p_dist,
       width = 10, height = 5, units = 'in', dpi = 600, bg = 'white')
