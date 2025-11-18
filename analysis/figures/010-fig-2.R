library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots with multiple panels
source('analysis/figures/000-default-ggplot-theme.R')

d <- tidyr::expand_grid(year = seq(1981, 2025, length.out = 400),
                        latitude = seq(-90, 90)) %>%
  mutate(mu_hat = rbeta(n(), 2, 2) * 1.1 - 0.1,
         denvar = runif(n(), 0, 0.04))

fig_2 <- plot_grid(
  ggplot() +
    geom_raster(aes(year, latitude, fill = mu_hat), d) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    scale_fill_gradientn('Draft mean NDVI', colours = ndvi_pal,
                         limits = c(-1, 1)) +
    theme(legend.position = 'top', legend.key.width = rel(2)),
  ggplot() +
    geom_raster(aes(year, latitude, fill = denvar), d) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_davos(name = 'Fake DENVar', limits = c(0, NA), reverse = TRUE) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    theme(legend.position = 'top', legend.key.width = rel(2)),
  labels = 'AUTO', ncol = 1)

ggsave('figures/fig-2.png', fig_2, width = 8.5, height = 10, dpi = 600,
       bg = 'white')
