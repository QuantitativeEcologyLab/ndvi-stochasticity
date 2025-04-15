library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for figures
library('khroma')  # for colorblind-friendly color scale
library('cowplot') # for ggplot plots in grids
source('analysis/figures/000-default-ggplot-theme.R')

r_mu <- rast('H:/GitHub/stochasticity_parks/Outputs/mean_predNDVI_raster.nc')
r_s2 <- rast('H:/GitHub/stochasticity_parks/Outputs/var_residuals_raster.nc')
r_dhi <- rast('data/other-rasters/dhi-data/dhi_ndvi_2015.tif')
r_hfi <- rast('data/other-rasters/hfi-layers/ml_hfi_v1_2019.nc')

layout(t(1:2))
plot(r_mu)
plot(r_s2)
layout(1)
plot(r_dhi)
plot(r_hfi)

d <- as.data.frame(r_s2, xy = TRUE) %>%
  rename(s2 = 3) %>%
  mutate(.,
         mu = extract(r_mu, select(., x, y), xy = FALSE)[, 2],
         hfi = extract(r_hfi, select(., x, y))[, 2]) %>%
  bind_cols(.,
            extract(r_dhi, select(., x, y)) %>%
              select(! ID) %>%
              rename(dhi_cumulative = dhi_ndvi_2015_1,
                     dhi_min = dhi_ndvi_2015_2,
                     dhi_seasonal = dhi_ndvi_2015_3)) %>%
  as_tibble()
d

YLIMS <- c(NA, ceiling(max(d$s2) * 100) / 100)
ZLIMS <- log(c(1, 150))

plot_hex <- function(variable, bins = 50, xlab, .data = d) {
  names(.data)[colnames(.data) == variable] <- 'variable'
  
  #' not adding `geom_smooth()` because it's too slow with large data sets
  .data %>%
    filter(, ! is.na(variable)) %>%
    ggplot() +
    geom_hex(aes(variable, s2, fill = log10(after_stat(count))),
             color = 'black', bins = bins, linewidth = 0.25)+
    scale_fill_acton(name = 'Count (log scale)', range = c(0, 1),
                     limits = ZLIMS, reverse = TRUE,
                     labels = \(.lab) round(exp(.lab))) +
    labs(x = xlab, y = 'DENVar') +
    ylim(YLIMS) +
    theme(legend.position = 'none') +
    #' for a safety check for fill limits, `ZLIMS`
    geom_hex(aes(variable, s2, color = log10(after_stat(count)) > ZLIMS[2],
                 lwd = log10(after_stat(count)) > ZLIMS[2]),
             fill = 'transparent', bins = bins, show.legend = FALSE) +
    scale_color_manual(values = c('#FFFFFF00', 'red'),
                       labels = c(FALSE, TRUE)) +
    scale_linewidth_manual(values = c(0.5, 10), labels = c(FALSE, TRUE))
}

plot_grid(
  plot_hex('mu', xlab = 'Mean NDVI') +
    theme(legend.position = 'top') + ggtitle('Full color vision'),
  colorblindr::cvd_grid(
    plot_hex('mu', xlab = 'Mean NDVI') +
      theme(legend.position = 'top')), ncol = 2)

p <- plot_grid(
  plot_hex('mu', xlab = 'Mean NDVI'),
  plot_hex('hfi', xlab = 'Human Footprint Index (ml-HFI)'),
  plot_hex('dhi_cumulative', xlab = 'Cumulative DHI'),
  plot_hex('dhi_min', xlab = 'Minimum DHI'),
  plot_hex('dhi_seasonal', xlab = 'Seasonal range in DHI'))

ggsave('figures/hexplots-canada-test.png', plot = p,
       width = 16, height = 8, units = 'in', dpi = 600, bg = 'white')
