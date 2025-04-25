library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots in grids
source('analysis/figures/000-default-ggplot-theme.R')

biomes <- rast('data/ecoregions/wwf-ecoregions.tif')
plot(biomes)
biomes_sf <- st_read('data/ecoregions/ecoregions-polygons.shp')
unique(biomes_sf$WWF_MHTNAM)
mu <- rast('data/output/mean-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif')
denvar <- rast('data/output/var-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif')
plot(denvar)

d <- as.data.frame(mu, xy = TRUE) %>%
  mutate(.,
         s2_hat = as.data.frame(denvar)[[1]],
         i = extract(biomes, select(., x, y))[, 2],
         biome = sort(unique(ecoregions$WWF_MHTNAM))[i + 1]) %>%
  filter(biome != 'Inland Water') %>%
  mutate(biome = factor(biome, levels = c(
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
  as_tibble()

biomes_sf <- biomes_sf %>%
  filter(WWF_MHTNAM != 'Inland Water') %>%
  mutate(WWF_MHTNAM = factor(WWF_MHTNAM, levels = levels(d$biome)))

biome_pal <- color('discreterainbow')(n_distinct(d$biome))

ggplot(d, aes(x, y, fill = biome)) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_discreterainbow(name = 'Biome') +
  theme(legend.position = 'top')

ggplot(d) +
  facet_wrap(~ biome, scales = 'free') +
  geom_histogram(aes(s2_hat, fill = biome), binwidth = 0.0005,
                 center = 0.0001, color = 'black') +
  labs(x = 'DENVar', y = 'Count') +
  scale_fill_discreterainbow() +
  theme(legend.position = 'none', strip.placement = 'outside',
        strip.text = element_text(size = rel(1)),
        strip.background = element_blank())

p_eco <-
  ggplot(biomes_sf) +
  geom_sf(aes(fill = WWF_MHTNAM), color = 'black', lwd = .05) +
  geom_hline(yintercept = 0, color = 'black', lwd = 0.1, lty = 'dashed') +
  scale_fill_discreterainbow(name = 'Ecoregion') +
  scale_x_continuous(expand = c(0, 0)) +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5),
        plot.margin = margin(5, 7.5, 5, 5))

make_hist <- function(.biome, var_only = FALSE) {
  .fill <- biome_pal[which(levels(d$biome) == .biome)]
  
  .title <- case_when(
    .biome == 'Tropical and Subtropical Grasslands, Savannas and Shrublands' ~
      'Tropical and Subtropical Grasslands,
Savannas and Shrublands',
    TRUE ~ .biome)
  
  sci_notation <- function(value) {
    .exp <- floor(log10(value))
    
    .labs <- paste0(floor(value / (10^(.exp))), 'e', .exp)
    
    if_else(value == 0, '0', .labs)
  }
  
  small_theme <-
    theme(legend.position = 'none', strip.placement = 'outside',
          strip.text = element_text(size = rel(0.5)),
          strip.background = element_blank(),
          plot.title = element_text(size = rel(0.5)),
          axis.text = element_text(size = rel(0.5)),
          axis.title = element_text(size = rel(0.5)))
  
  if(var_only) {
    hists <-
      d %>%
      filter(biome == .biome) %>%
      ggplot() +
      geom_histogram(aes(s2_hat), binwidth = 0.002, center = 0.0001,
                     fill = .fill, color = 'transparent', na.rm = TRUE) +
      labs(x = 'DENVar', y = NULL) +
      scale_x_continuous(limits = range(d$s2_hat), breaks = c(0, 0.05, 0.1)) +
      scale_y_continuous(labels = sci_notation) +
      scale_fill_discreterainbow(breaks = levels(d$biome)) +
      small_theme
  } else {
    hists <-
      plot_grid(
        d %>%
          filter(biome == .biome) %>%
          ggplot() +
          geom_histogram(aes(mu_hat), binwidth = 0.02, center = 0.01,
                         fill = .fill, color = 'transparent', na.rm = TRUE) +
          labs(x = 'Mean NDVI', y = NULL) +
          scale_x_continuous(limits = c(-0.1, 0.45), breaks = c(0, 0.2, 0.4)) +
          scale_y_continuous(labels = sci_notation) +
          scale_fill_discreterainbow(breaks = levels(d$biome)) +
          small_theme,
        d %>%
          filter(biome == .biome) %>%
          ggplot() +
          geom_histogram(aes(s2_hat), binwidth = 0.002, center = 0.0001,
                         fill = .fill, color = 'transparent', na.rm = TRUE) +
          labs(x = 'DENVar', y = NULL) +
          scale_x_continuous(limits = range(d$s2_hat), breaks = c(0, 0.05, 0.1)) +
          scale_y_continuous(labels = sci_notation) +
          scale_fill_discreterainbow(breaks = levels(d$biome)) +
          small_theme)
  }
  
  (plot_grid(
    ggdraw() +
      draw_label(.title, fontface = 'bold', x = 0.5,
                 hjust = 0.5, size = 5, vjust = 0.3),
    hists, ncol = 1, rel_heights = c(1, 10)) +
      theme(plot.background = element_rect(color = 'grey', linewidth = 0.5),
            plot.margin = margin(5, 5, 5, 5))) %>%
    plot_grid(NULL, ., NULL, nrow = 1, rel_widths = c(1, 40, 1)) %>%
    plot_grid(NULL, ., NULL, ncol = 1, rel_heights = c(1, 20, 1))
}

make_hist(.biome = 'Tropical and Subtropical Grasslands, Savannas and Shrublands')
make_hist(.biome = 'Tropical and Subtropical Grasslands, Savannas and Shrublands',
          var_only = TRUE)

p <-
  plot_grid(
    plot_grid(
      plot_grid(plotlist = map(levels(d$biome)[1:4], make_hist), ncol = 1),
      p_eco,
      plot_grid(plotlist = map(levels(d$biome)[12:15], make_hist), ncol = 1),
      rel_widths = c(1, 5, 1), nrow = 1),
    plot_grid(plotlist = map(levels(d$biome)[5:11], make_hist), nrow = 1),
    ncol = 1, rel_heights = c(4, 1))

ggsave('figures/global-models/ecoregions-histograms.png', p,
       width = 10.1 * 7/5, height = 6 * 5/4, units = 'in',
       dpi = 600, bg = 'white')

p <-
  plot_grid(
    plot_grid(
      plot_grid(plotlist = map(levels(d$biome)[1:4],
                               \(.b) make_hist(.b, TRUE)), ncol = 1),
      p_eco,
      plot_grid(plotlist = map(levels(d$biome)[12:15],
                               \(.b) make_hist(.b, TRUE)), ncol = 1),
      rel_widths = c(1, 5, 1), nrow = 1),
    plot_grid(plotlist = map(levels(d$biome)[5:11],
                             \(.b) make_hist(.b, TRUE)), nrow = 1),
    ncol = 1, rel_heights = c(4, 1))

ggsave('figures/global-models/ecoregions-histograms-var-only.png', p,
       width = 10.1 * 7/5, height = 6 * 5/4, units = 'in',
       dpi = 600, bg = 'white')
