library('dplyr')     # for data wrangling
library('sf')        # for simple features
library('terra')     # for rasters
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
library('tidyterra') # for plotting spatRasters with ggplot2
library('khroma')    # for colorblind-friendly palettes
source('analysis/figures/000-default-ggplot-theme.R')
color('okabeito')(8)[c(3, 4, 6, 7)]
pal_groups <- c('#56B4E9', 'grey', '#009E73', '#0072B2', '#D55E00')
plot_scheme_colorblind(pal_groups[c(1, 2, 4, 3, 5)]) # achromatic view ok

robinson_crs <- '+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

shp <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  mutate(WWF_MHTNAM = factor(WWF_MHTNAM, levels = c(
    'Rock and Ice',
    'Tundra',
    'Boreal Forests/Taiga',
    'Montane Grasslands and Shrublands',
    'Temperate Grasslands, Savannas and Shrublands',
    'Temperate Broadleaf and Mixed Forests',
    'Temperate Conifer Forests',
    'Mediterranean Forests, Woodlands and Scrub',
    'Mangroves',
    'Tropical and Subtropical Moist Broadleaf Forests',
    'Tropical and Subtropical Dry Broadleaf Forests',
    'Tropical and Subtropical Coniferous Forests',
    'Tropical and Subtropical Grasslands, Savannas and Shrublands',
    'Flooded Grasslands and Savannas',
    'Deserts and Xeric Shrublands'))) %>%
  st_transform(robinson_crs)

bounds <- tibble(x = c(-180:180, rep(180, 181), 180:-180, rep(-180, 181)),
                 y = c(rep(-90, 361), -90:90, rep(90, 361), 90:-90)) %>%
  st_as_sf(coords = c('x', 'y'), crs = 'EPSG:4326') %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast('POLYGON') %>%
  st_transform(robinson_crs)

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
  
  plot_scheme(color('devon')(10))
  plot_scheme(color('lajolla')(10)) # for elevation
  plot_scheme(color('bamako')(10)) # for n rasters
  plot_scheme(color('davos')(10)) # for DENVar
  plot_scheme(color('bilbao')(10))
  plot_scheme(color('nuuk')(10))
  plot_scheme(color('oslo')(10))
  plot_scheme(color('grayC')(10))
  plot_scheme(color('hawaii')(10))
  plot_scheme(color('lapaz')(10))
  plot_scheme(color('tokyo')(10))
  plot_scheme(color('buda')(10))
  plot_scheme(color('acton')(10))
  plot_scheme(color('turku')(10))
  plot_scheme(color('imola')(10))
  plot_scheme(color('batlow')(10))
  plot_scheme(color('batlowW')(10))
  plot_scheme(color('batlowK')(10))
  plot_scheme(color('brocO')(10))
  plot_scheme(color('corkO')(10))
  plot_scheme(color('vikO')(10))
  plot_scheme(color('romaO')(10))
  plot_scheme(color('bamO')(10))
  plot_scheme(color('YlOrBr')(10))
  plot_scheme(color('iridescent')(10))
  plot_scheme(color('incandescent')(10))
  plot_scheme(color('smoothrainbow')(10))
}

# biomes ----
n_distinct(shp$WWF_MHTNAM)

p_biome <-
  ggplot(shp) +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(aes(fill = WWF_MHTNAM), color = 'black', lwd = .05) +
  geom_hline(yintercept = 0, color = 'black', lwd = 0.1, lty = 'dashed') +
  scale_fill_discreterainbow(name = 'Biome') +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5),
        panel.grid = element_blank(), text = element_text(face = 'bold'))

ggsave('figures/input-data/biomes.png', p_biome,
       width = 10, height = 6.5, units = 'in', dpi = 1e3, bg = 'white')

# data groups ----
n_distinct(shp$group)

p_group <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(aes(fill = group), shp, color = 'black', lwd = .05) +
  geom_hline(yintercept = 0, color = 'black', lwd = 0.1, lty = 'dashed') +
  scale_fill_manual('Group', values = pal_groups) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top', legend.text = element_text(size = 8),
        panel.grid = element_blank(), text = element_text(face = 'bold'))

ggsave('figures/input-data/data-groups.png', p_group,
       width = 10, height = 6, units = 'in', dpi = 1e3, bg = 'white')

# elevation ----
r_elev <- rast('data/elev-raster.tif') %>%
  crop(., st_transform(shp, crs(.)), mask = TRUE, touches = FALSE) %>%
  project(robinson_crs)

# qattara depression and dead sea are < 0 m but also separate polygons
ggplot() +
  geom_spatraster(data = r_elev < 0) +
  scale_fill_manual(values = c('white', 'blue'))

# the only elevations below minimum dead sea elevation are coastlines
ggplot() +
  geom_spatraster(data = r_elev < -432) +
  scale_fill_manual(values = c('white', 'red'))

hist(r_elev, breaks = 1000)
abline(v = -432, col = 'red')
mean(values(r_elev < -432), na.rm = TRUE)

# change coastlines to 0 for the plot
values(r_elev) <- if_else(values(r_elev) < -432, -432, values(r_elev))

p_elev <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_spatraster(data = r_elev) +
  geom_sf(data = shp, color = 'black', fill = 'transparent', lwd = 0.05) +
  scale_fill_lajolla(name = 'Elevation above sea level (m)') +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5),
        panel.grid = element_blank(), text = element_text(face = 'bold'))

ggsave('figures/input-data/elev-m.png', p_elev,
       width = 10, height = 6, units = 'in', dpi = 1e3, bg = 'white')

#' function similar to `khroma::plot_scheme_colorblind()`
plot_pal <- function(pal) {
  names(pal) <- NULL
  
  .p <- colorblindr::palette_plot(pal, color_labels = FALSE,
                                  xmargin = 0)
  
  plot_grid(
    # figure in full trichromatic vision
    .p +
      ggtitle('Trichromatic vision') +
      theme(plot.title = element_text(face = 'bold')),
    # figure in three color-vision deficient visions
    colorblindr::cvd_grid(.p, severity = 1), ncol = 1,
    rel_heights = c(1.25, 2))
}

plot_pal(pal_groups)

p_palettes <- plot_grid(
  NULL, NULL, NULL,
  plot_pal(ndvi_pal),
  plot_pal(rev(color('davos')(1e3))),
  plot_pal(rev(color('bamako')(3))), # A2
  NULL, NULL, NULL,
  plot_pal(pal_groups), # A3
  plot_pal(color('discreterainbow')(n_distinct(shp$WWF_MHTNAM))), # A4
  plot_pal(color('lajolla')(1e3)), # A5
  label_x = 0.5, hjust = 0.5, nrow = 4, rel_heights = c(0.03, 1),
  labels = c('NDVI', 'DENVar', 'Fig. A2', '', '', '', 
             paste0('Fig. A', 3:5), '', '', ''))

ggsave('figures/input-data/color-palettes.png', p_palettes,
       width = 8, height = 6, scale = 2, units = 'in', dpi = 300,
       bg = 'white')
