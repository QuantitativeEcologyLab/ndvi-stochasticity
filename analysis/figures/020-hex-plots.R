library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for figures
library('khroma')  # for colorblind-friendly color scale
library('cowplot') # for ggplot plots in grids
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')

m_mu <- readRDS('models/global-test/test-mean-gam-sos-only.rds')
m_s2 <- readRDS('models/global-test/test-var-gam-sos-only.rds')

r_dhi <- rast('data/other-rasters/dhi-data/dhi_ndvi_2015.tif')
r_hfi <- rast('data/other-rasters/hfi-layers/ml_hfi_v1_2019.nc')
r_rich <- rast('data/other-rasters/iucn-red-list-spp-richness/Combined_SR_2024.tif') %>%
  project('EPSG:4326') # all others are lat/long
crs(r_dhi, proj = TRUE) == crs(r_hfi, proj = TRUE) &
  crs(r_hfi, proj = TRUE) == crs(r_rich, proj = TRUE)

# check rasters
plot(r_dhi)
plot(r_hfi)
plot(r_rich)

#' *will need to change the data later*
d <- m_mu$model %>%
  mutate(mu_hat = fitted(m_mu),
         s2_hat = fitted(m_s2),
         hfi = extract(r_hfi, select(., x, y))[, 2],
         richness = extract(r_rich, select(., x, y))[, 2]) %>%
  bind_cols(.,
            extract(r_dhi, select(., x, y)) %>%
              select(! ID) %>%
              rename(dhi_cumulative = dhi_ndvi_2015_1,
                     dhi_min = dhi_ndvi_2015_2,
                     dhi_seasonal = dhi_ndvi_2015_3)) %>%
  as_tibble()
d

d %>%
  tidyr::pivot_longer(c(mu_hat, hfi:dhi_seasonal), names_to = 'variable',
                      values_to = 'value') %>%
  mutate(variable = case_when(variable == 'mu_hat' ~ 'Mean NDVI',
                              variable == 'hfi' ~ 'ml-HFI',
                              variable == 'richness' ~ 'Species richness',
                              variable == 'dhi_cumulative' ~ 'Cumulative DHI',
                              variable == 'dhi_min' ~ 'Minimum DHI',
                              variable == 'dhi_seasonal' ~ 'Seasonal range in DHI') %>%
           factor(., levels = unique(.))) %>%
  ggplot() +
  facet_wrap(~ variable, scales = 'free_x', strip.position = 'bottom') +
  geom_hex(aes(value, s2_hat, fill = log10(after_stat(count))),
           color = 'black', bins = 20, linewidth = 0.1, na.rm = TRUE) +
  scale_fill_iridescent(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, breaks = seq(0, 6, by = 2),
    labels = \(.x) paste0('1e', .x)) +
  labs(x = NULL, y = 'DENVar') +
  theme(legend.position = 'top', strip.background = element_blank(),
        strip.placement = 'outside',
        strip.text = element_text(size = rel(1)))

ggsave('figures/hexplots-canada-test.png',
       width = 12, height = 10, units = 'in', dpi = 600, bg = 'white')

#' YLIMS <- c(NA, ceiling(max(d$s2_hat) * 100) / 100)
#' ZLIMS <- log10(c(1, 3e4))
#' 
#' plot_hex <- function(variable, bins = 100, xlab = 'x', .data = d) {
#'   names(.data)[colnames(.data) == variable] <- 'variable'
#'   
#'   #' not adding `geom_smooth()` because it's too slow with large data sets
#'   .data %>%
#'     filter(, ! is.na(variable)) %>%
#'     ggplot() +
#'     geom_hex(aes(variable, s2_hat, fill = log10(after_stat(count))),
#'              color = 'black', bins = bins, linewidth = 0.1) +
#'     scale_fill_iridescent(
#'       name = expression(paste(bold('Count (log'), bold(''['10']),
#'                               bold(' scale)'))), range = c(0, 1),
#'       limits = ZLIMS, reverse = FALSE, breaks = seq(0, 6, by = 2),
#'       labels = \(.x) paste0('1e', .x)) +
#'     labs(x = xlab, y = 'DENVar') +
#'     ylim(YLIMS) +
#'     theme(legend.position = 'none') +
#'     #' for a safety check for fill limits, `ZLIMS`
#'     geom_hex(aes(variable, s2_hat,
#'                  color = log10(after_stat(count)) > ZLIMS[2],
#'                  lwd = log10(after_stat(count)) > ZLIMS[2]),
#'              fill = 'transparent', bins = bins, show.legend = FALSE) +
#'     scale_color_manual(values = c('transparent', 'red'),
#'                        labels = c(FALSE, TRUE)) +
#'     scale_linewidth_manual(values = c(0.5, 3), labels = c(FALSE, TRUE))
#' }
#' 
#' plot_hex('mu_hat', xlab = 'Mean NDVI') +
#'   theme(legend.position = 'top')
#' 
#' khroma::info()
#' 
#' plot_grid(
#'   plot_hex('mu_hat', xlab = 'Mean NDVI') +
#'     theme(legend.position = 'top') +
#'     ggtitle('Full color vision'),
#'   colorblindr::cvd_grid(plot_hex('mu_hat', xlab = 'Mean NDVI') +
#'                           theme(legend.position = 'top')),
#'   ncol = 2)
#' 
#' p <- plot_grid(
#'   plot_hex('mu_hat', xlab = 'Mean NDVI'),
#'   plot_hex('hfi', xlab = 'Human Footprint Index (ml-HFI)'),
#'   plot_hex('richness', xlab = 'Species richness'),
#'   plot_hex('dhi_cumulative', xlab = 'Cumulative DHI'),
#'   plot_hex('dhi_min', xlab = 'Minimum DHI'),
#'   plot_hex('dhi_seasonal', xlab = 'Seasonal range in DHI')) %>%
#'   plot_grid(get_legend(plot_hex('mu_hat') +
#'                          theme(legend.position = 'top',
#'                                legend.key.width = unit(1, 'in'))),
#'             ., rel_heights = c(1, 10), ncol = 1, labels = 'AUTO')
#' p
#' 
#' ggsave('figures/hexplots-canada-test.png', plot = p,
#'        width = 12, height = 6, units = 'in', dpi = 600, bg = 'white')
