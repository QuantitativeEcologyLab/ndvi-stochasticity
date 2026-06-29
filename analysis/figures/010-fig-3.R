library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots with multiple panels
source('analysis/figures/000-default-ggplot-theme.R')

d <- readr::read_csv(d, 'output/yearly-estimates.csv', num_threads = 10)

d <- d %>%
  summarize(mu_hat = mean(mu_hat, na.rm = TRUE),
            s2_hat = mean(s2_hat, na.rm = TRUE),
            .by = c(year, y))

mean(d$s2_hat > 0.05)
d <- mutate(d,
            s2_hat = if_else(s2_hat > 0.05, 0.05, s2_hat),
            ref_mu = mean(mu_hat[year <= 1990]),
            ref_s2 = mean(s2_hat[year <= 1990]),
            diff_mu = mu_hat - ref_mu,
            diff_s2 = s2_hat - ref_s2,
            .by = y)

fig_3 <- plot_grid(
  ggplot(d, aes(year, y)) +
    geom_raster(aes(fill = mu_hat)) +
    geom_contour(aes(year, y, z = mu_hat), color = '#00000080', bins = 8) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    scale_fill_gradientn('Mean NDVI   ', colours = ndvi_pal,
                         limits = c(-1, 1)) +
    theme(legend.position = 'top', legend.key.width = rel(2),
          legend.title = element_text(hjust = 0)),
  ggplot(d, aes(year, y)) +
    geom_raster(aes(fill = s2_hat)) +
    geom_contour(aes(z = s2_hat), color = '#00000080', bins = 8) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_davos(name = 'DENVar ', limits = c(0, NA), reverse = TRUE) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    theme(legend.position = 'top', legend.key.width = rel(2),
          legend.title = element_text(hjust = 0)),
  ggplot(d, aes(year, y)) +
    geom_raster(aes(fill = diff_mu)) +
   geom_contour(aes(year, y, z = diff_mu), color = '#00000080', bins = 8) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    scale_fill_vik(name = 'Change in mean NDVI      ', reverse = TRUE) +
    theme(legend.position = 'top', legend.key.width = rel(2),
          legend.title = element_text(hjust = 0)),
  ggplot(d, aes(year, y)) +
    geom_raster(aes(fill = diff_s2)) +
   geom_contour(aes(z = diff_s2), color = '#00000080', bins = 8) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_bam(name = 'Change in DENVar    ', reverse = TRUE) +
    labs(x = 'Year', y = 'Latitude (degrees N)') +
    theme(legend.position = 'top', legend.key.width = rel(2),
          legend.title = element_text(hjust = 0)),
  labels = 'AUTO', ncol = 2)

ggsave('figures/fig-3.png', fig_3, width = 17, height = 10, dpi = 600,
       bg = 'white')
