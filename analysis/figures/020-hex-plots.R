library('sf')      # for simple features
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for figures
library('khroma')  # for colorblind-friendly color scale
library('cowplot') # for ggplot plots in grids
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')

r_s2 <- rast('data/output/var-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif')
ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  # drop polygons not in data
  filter(poly_id %in% names(readRDS('data/ecoregions/poly-nbs-global.rds'))) %>%
  vect() %>% # convert to spatVect
  project(crs(r_s2))
r_s2 <- crop(r_s2, ecoregions, mask = TRUE)
r_mu <- rast('data/output/mean-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif')
r_dhi <- rast('data/other-rasters/dhi-data/dhi_ndvi_2015.tif')
r_hfi <- rast('data/other-rasters/hfi-layers/ml_hfi_v1_2019.nc') # no data on greenland or arctic
r_rich <- rast('data/other-rasters/iucn-red-list-spp-richness/Combined_SR_2024.tif') %>%
  project(crs(r_s2))
n_distinct(sapply(list(r_s2, r_mu, r_dhi, r_hfi, r_rich), \(x) crs(x, proj = TRUE)))

# check rasters
if(FALSE) {
  plot(r_s2)
  plot(r_mu)
  plot(r_dhi)
  plot(r_hfi)
  plot(r_rich)
}

d <- as.data.frame(r_mu, xy = TRUE) %>%
  mutate(.,
         s2_hat = extract(r_s2, select(., x, y))[, 2],
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
  filter(is.na(hfi)) %>%
  plot(y ~ x, .)

d %>%
  filter(mu_hat > -0.25) %>% #' *remove: there were 4 points mean < -0.25*
  tidyr::pivot_longer(c(mu_hat, hfi:dhi_seasonal), names_to = 'variable',
                      values_to = 'value') %>%
  mutate(lab = case_when(variable == 'mu_hat' ~ 'Mean NDVI',
                         variable == 'hfi' ~ 'ml-HFI',
                         variable == 'richness' ~ 'Species richness',
                         variable == 'dhi_cumulative' ~ 'Cumulative DHI',
                         variable == 'dhi_min' ~ 'Minimum DHI',
                         variable == 'dhi_seasonal' ~ 'Seasonal range in DHI') %>%
           factor(., levels = unique(.))) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', strip.position = 'left') +
  geom_hex(aes(s2_hat, value, fill = log10(after_stat(count))),
           color = 'black', bins = 50, linewidth = 0.1, na.rm = TRUE) +
  scale_fill_iridescent(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, labels = \(.x) 10^.x) +
  labs(y = NULL, x = 'DENVar') +
  theme(legend.position = 'top', strip.background = element_blank(),
        strip.placement = 'outside', legend.key.width = unit(0.5, 'in'),
        strip.text = element_text(size = rel(1)))

ggsave('figures/hexplots-global-test.png',
       width = 16, height = 10, units = 'in', dpi = 600, bg = 'white')
