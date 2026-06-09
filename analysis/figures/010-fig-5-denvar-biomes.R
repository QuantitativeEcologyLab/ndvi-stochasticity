library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots in grids
source('analysis/figures/000-default-ggplot-theme.R')
source('analysis/figures/000-robinson-objects.R')

d <- readRDS('output/d_rob.rds')

biomes_sf <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  st_transform(robinson_crs)
biomes <- rasterize(biomes_sf, rast(d), field = 'biome', touches = TRUE)

plot(biomes)

d <- d %>%
  mutate(biome = factor(biome, levels = c(
    'Rock and Ice',
    'Tundra',
    'Boreal Forests/Taiga',
    'Montane Grasslands and Shrublands',
    'Temperate Grasslands, Savannas and Shrublands',
    'Temperate Broadleaf and Mixed Forests',
    'Temperate Coniferous Forests',
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
  filter(biome != 'Inland Water') %>%
  mutate(biome = factor(biome, levels = levels(d$biome)))

biome_pal <- color('discreterainbow')(n_distinct(d$biome))

ggplot(d) +
  facet_wrap(~ biome, scales = 'free') +
  geom_histogram(aes(s2_hat, fill = biome), binwidth = 0.0005,
                 center = 0.0001, color = 'black') +
  labs(x = 'DENVar', y = 'Count') +
  scale_fill_discreterainbow() +
  theme(legend.position = 'none', strip.placement = 'outside',
        strip.text = element_text(size = rel(1)),
        strip.background = element_blank())

p_biome <-
  ggplot(biomes_sf) +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(aes(fill = biome), color = 'black', lwd = .05) +
  geom_hline(yintercept = 0, color = 'black', lwd = 0.1, lty = 'dashed') +
  scale_fill_discreterainbow(name = 'Biome') +
  scale_x_continuous(expand = c(0, 0)) +
  theme_void() +
  theme(legend.position = 'top', legend.text = element_text(size = 4.5),
        panel.grid = element_blank(), text = element_text(face = 'bold'))

marginal_plot <- function(.biome, var_only = FALSE, hex = TRUE) {
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
  
  if(hex) {
    if(var_only) {
      warning('`var_only` only works when `hex` is FALSE.')
    }
    p <-
      d %>%
      filter(biome == .biome) %>%
      ggplot() +
      geom_hex(aes(mu_hat, s2_hat), na.rm = TRUE) +
      labs(x = 'Mean NDVI', y = 'DENVar') +
      scale_y_continuous(limits = range(d$s2_hat),
                         breaks = c(0, 0.02, 0.04, 0.06)) +
      scale_fill_lapaz(reverse = TRUE) +
      small_theme
  } else {
    if(var_only) {
      p <-
        d %>%
        filter(biome == .biome) %>%
        ggplot() +
        geom_histogram(aes(s2_hat), binwidth = 0.002, center = 0.0001,
                       fill = .fill, color = 'transparent', na.rm = TRUE) +
        labs(x = 'DENVar', y = NULL) +
        scale_x_continuous(limits = range(d$s2_hat),
                           breaks = c(0, 0.02, 0.04, 0.06)) +
        scale_y_continuous(labels = sci_notation) +
        scale_fill_discreterainbow(breaks = levels(d$biome)) +
        small_theme
    } else {
      p <-
        plot_grid(
          d %>%
            filter(biome == .biome) %>%
            ggplot() +
            geom_histogram(aes(mu_hat), binwidth = 0.02, center = 0.01,
                           fill = .fill, color = 'transparent', na.rm = TRUE) +
            labs(x = 'Mean NDVI', y = NULL) +
            scale_x_continuous(limits = c(-0.1, max(d$mu_hat, na.rm = TRUE)),
                               breaks = seq(0, 1, by = 0.2)) +
            scale_y_continuous(labels = sci_notation) +
            scale_fill_discreterainbow(breaks = levels(d$biome)) +
            small_theme,
          d %>%
            filter(biome == .biome) %>%
            ggplot() +
            geom_histogram(aes(s2_hat), binwidth = 0.002, center = 0.0001,
                           fill = .fill, color = 'transparent', na.rm = TRUE) +
            labs(x = 'DENVar', y = NULL) +
            scale_x_continuous(limits = range(d$s2_hat),
                               breaks = c(0, 0.02, 0.04, 0.06)) +
            scale_y_continuous(labels = sci_notation) +
            scale_fill_discreterainbow(breaks = levels(d$biome)) +
            small_theme)
    }
  }
  
  p <- plot_grid(
    ggdraw() +
      draw_label(.title, fontface = 'bold', x = 0.2,
                 hjust = -0.2, size = 5, vjust = 0.3),
    p, ncol = 1, rel_heights = c(1, 10))
  
  if(! hex & ! var_only) {
    p <- p +
      theme(plot.background = element_rect(color = 'grey', linewidth = 0.5),
            plot.margin = margin(5, 5, 5, 5))
  }
  
  plot_grid(NULL, p, NULL, nrow = 1, rel_widths = c(1, 40, 1)) %>%
    plot_grid(NULL, ., NULL, ncol = 1, rel_heights = c(1, 20, 1))
}

marginal_plot(.biome = 'Tropical and Subtropical Grasslands, Savannas and Shrublands')
marginal_plot(.biome = 'Tropical and Subtropical Grasslands, Savannas and Shrublands',
              hex = FALSE)
marginal_plot(.biome = 'Tropical and Subtropical Grasslands, Savannas and Shrublands',
              hex = FALSE, var_only = TRUE)

p_hex <-
  plot_grid(
    plot_grid(
      plot_grid(plotlist = map(levels(d$biome)[1:4], marginal_plot), ncol = 1),
      p_biome,
      plot_grid(plotlist = map(levels(d$biome)[12:15], marginal_plot), ncol = 1),
      rel_widths = c(1, 5, 1), nrow = 1),
    plot_grid(plotlist = map(levels(d$biome)[5:11], marginal_plot), nrow = 1),
    ncol = 1, rel_heights = c(4, 1))

ggsave('figures/fig-5.png', p_hex,
       width = 10.1 * 7/5, height = 6 * 5/4, units = 'in',
       dpi = 600, bg = 'white')

p_hist <-
  plot_grid(
    plot_grid(
      plot_grid(plotlist = map(levels(d$biome)[1:4],
                               \(.b) marginal_plot(.b, FALSE, hex = FALSE)),
                ncol = 1),
      p_biome,
      plot_grid(plotlist = map(levels(d$biome)[12:15],
                               \(.b) marginal_plot(.b, FALSE, hex = FALSE)),
                ncol = 1),
      rel_widths = c(1, 5, 1), nrow = 1),
    plot_grid(plotlist = map(levels(d$biome)[5:11],
                             \(.b) marginal_plot(.b, FALSE, hex = FALSE)),
              nrow = 1),
    ncol = 1, rel_heights = c(4, 1))

ggsave('figures/global-models/biomes-histograms.png', p_hist,
       width = 10.1 * 7/5, height = 6 * 5/4, units = 'in',
       dpi = 600, bg = 'white')

p_hist_var <-
  plot_grid(
    plot_grid(
      plot_grid(plotlist = map(levels(d$biome)[1:4],
                               \(.b) marginal_plot(.b, TRUE, hex = FALSE)),
                ncol = 1),
      p_biome,
      plot_grid(plotlist = map(levels(d$biome)[12:15],
                               \(.b) marginal_plot(.b, TRUE, hex = FALSE)),
                ncol = 1),
      rel_widths = c(1, 5, 1), nrow = 1),
    plot_grid(plotlist = map(levels(d$biome)[5:11],
                             \(.b) marginal_plot(.b, TRUE, hex = FALSE)),
              nrow = 1),
    ncol = 1, rel_heights = c(4, 1))

ggsave('figures/global-models/biomes-histograms-var-only.png', p_hist_var,
       width = 10.1 * 7/5, height = 6 * 5/4, units = 'in',
       dpi = 600, bg = 'white')
