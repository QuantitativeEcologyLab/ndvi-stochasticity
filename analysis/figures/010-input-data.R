library('dplyr')     # for data wrangling
library('sf')        # for simple features
library('terra')     # for rasters
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
library('tidyterra') # for plotting spatRasters with ggplot2
library('khroma')    # for colorblind-friendly palettes
source('analysis/figures/000-default-ggplot-theme.R')
pal_polys <- c(color('okabeito')(8), 'grey')

eco <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(WWF_MHTNAM != 'Inland Water') %>%
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
    'Deserts and Xeric Shrublands')))

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
  plot_scheme(color('lajolla')(10))
  plot_scheme(color('bamako')(10))
  plot_scheme(color('davos')(10))
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

# individual polygons ----
n_distinct(eco$poly_id)
plot_scheme_colorblind(pal_polys)
ggsave('polygons-palette.png', path = 'figures/input-data',
       width = 6, height = 4, units = 'in', dpi = 300, bg = 'white')

set.seed(-3)
p_poly <-
  ggplot(eco) +
  geom_sf(fill = sample(pal_polys, nrow(eco), replace = TRUE),
          lwd = .15, color = 'white') +
  geom_sf(fill = 'transparent', lwd = .05, color = 'black') +
  scale_x_continuous(expand = c(0, 0))

ggsave('figures/input-data/polygons.png', p_poly,
       width = 10, height = 4.4, units = 'in', dpi = 1200, bg = 'white')

# ecoregions ----
n_distinct(eco$WWF_MHTNAM)

p_biome <-
  ggplot(eco) +
  geom_sf(aes(fill = WWF_MHTNAM), color = 'black', lwd = .05) +
  geom_hline(yintercept = 0, color = 'black', lwd = 0.1, lty = 'dashed') +
  scale_fill_discreterainbow(name = 'Biome') +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5))

ggsave('figures/input-data/ecoregions.png', p_biome,
       width = 10.1, height = 6, units = 'in', dpi = 1e3, bg = 'white')

# elevation ----
r_elev <- rast('data/elev-raster.tif') %>%
  crop(., st_transform(eco, crs(.)), mask = TRUE)

# change coastal values to 0 m above sea level
ggplot() +
  geom_spatraster(data = r_elev < 0)

# qattara depression and dead sea are < 0 m but also separate polygons
ggplot() +
  geom_spatraster(data = r_elev < 0) +
  scale_fill_manual(values = c('white', 'red'))

hist(r_elev, breaks = 1000)
abline(v = 0, col = 'red')
mean(values(r_elev < 0), na.rm = TRUE)

values(r_elev) <- if_else(values(r_elev) < 0, 0, values(r_elev))

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

# distance from coast
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

# hex plot of rasters over day of year and year
dates <- tibble(
  file_name = list.files('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files',
                         pattern = '*.nc', recursive = FALSE),
  date = substr(file_name,
                nchar(file_name) - nchar('YYYYmmdd_xYYYYmmddHHMMSS.nc') + 1,
                nchar(file_name) - nchar('_xYYYYmmddHHMMSS.nc')) %>%
    as.Date(format = '%Y%m%d'),
  doy = lubridate::yday(date),
  year = lubridate::year(date))

p_n <-
  dates %>%
  ggplot(aes(year, doy)) +
  coord_equal(ratio = 1 / 15) +
  geom_bin_2d(binwidth = c(0.9999999, 15)) +
  xlab('Year') +
  scale_y_continuous('Day of year', expand = c(0, 0),
                     breaks = c(1, 100, 200, 300, 366)) +
  scale_fill_lapaz(name = 'Number of rasters', reverse = TRUE,
                   limits = c(1, NA), range = c(0, 1),
                   breaks = c(1, 5, 10, 15)) +
  theme(legend.position = 'top')
p_n

ggsave('figures/input-data/n-rasters-time.png', p_n,
       width = 10, height = 5, units = 'in', dpi = 300, bg = 'white')

#' function similar to `khroma::plot_scheme_colorblind()`
plot_pal <- function(pal) {
  .d <- tibble(x = floor((1:length(pal))),
               y = x %% 2)
  
  .p <- colorblindr::palette_plot(pal, color_labels = FALSE,
                                  xmargin = 0)
  
  plot_grid(.p +
              ggtitle('Trichromatic') +
              theme(plot.title = element_text(face = 'bold')),
            colorblindr::cvd_grid(.p, severity = 1), ncol = 1,
            rel_heights = c(1.25, 2))
}

plot_pal(pal_polys)

plot_grid(NULL,
          plot_pal(rev(color('lapaz')(1e3))),
          NULL,
          plot_pal(color('discreterainbow')(n_distinct(eco$WWF_MHTNAM))),
          NULL,
          plot_pal(pal_polys),
          NULL,
          plot_pal(color('lajolla')(1e3)),
          # NULL,
          # plot_pal(color('davos')(1e3)),
          label_x = 0, ncol = 4,  rel_widths = c(0.075, 1),
          labels = c('A2.', '', 'A3.', '', 'A4.', '', 'A5.'))

ggsave('figures/input-data/color-palettes.png', width = 8, height = 6,
       scale = 2, units = 'in', dpi = 300, bg = 'white')
