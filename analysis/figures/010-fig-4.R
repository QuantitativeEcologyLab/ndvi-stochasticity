library('sf')      # for simple features
library('tidyr')   # for pivoting and nesting data
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for figures
library('khroma')  # for colorblind-friendly color scale
library('cowplot') # for ggplot plots in grids
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')

r_mu <- rast('output/long-term-estimates.tif')[[1]]
r_s2 <- rast('output/long-term-estimates.tif')[[2]]
ecoregions <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  vect() %>% # convert to spatVect
  project(crs(r_s2))
r_dhi <- rast('data/other-rasters/dhi-data/dhi_ndvi_2015.tif') %>%
  project(r_s2)
# using Keys HFI raster because it's faster; other one has different CRS
r_hfi <- rast('data/other-rasters/hfi-layers/ml_hfi_v1_2019.nc') %>%
  project(r_s2)
# r_hfi <- rast('H:/GitHub/ndvi-stochasticity/data/other-rasters/hfp_2021_100m_v1-2_cog.tif')
r_rich <- rast('data/other-rasters/iucn-red-list-spp-richness/Combined_SR_2024.tif') %>%
  project(crs(r_s2))
if(file.exists('data/other-rasters/FireCCI51based/global_monthly_burned_area_fraction_05deg_caluclated-average.nc')) {
  r_burn <- rast('data/other-rasters/FireCCI51based/global_monthly_burned_area_fraction_05deg_caluclated-average.nc')
} else {
  r_burn <- list.files('data/other-rasters/FireCCI51based',
                       pattern = '.nc', full.names = TRUE) %>%
    lapply(rast) %>%
    rast() %>%
    mean() %>%
    project(crs(ecoregions)) %>%
    mask(ecoregions)
  writeCDF(r_burn, 'data/other-rasters/FireCCI51based/global_monthly_burned_area_fraction_05deg_caluclated-average.nc')
}

r_precip <- rast('data/precipitation-yearly-mean-m-day.tif') %>%
  project(r_s2) %>%
  mask(ecoregions)

r_dcv <- rast('https://github.com/ahyoung-lim/Arbo_riskmaps_public/raw/refs/heads/main/outputs/Rasters/DCZ_riskmap_wmean_masked.tif')

#' *data that might be worth adding:*
#' - temperature (need to make a raster): `https://cds.climate.copernicus.eu/datasets/derived-near-surface-meteorological-variables?tab=download`
#' - seasonal temperature range (need to make a raster)
#' - max body size (can't find a raster)
#' - invasive species (can't find a raster)

# check rasters
if(FALSE) {
  plot(r_s2)
  plot(r_mu)
  plot(r_dhi)
  plot(r_hfi)
  plot(r_rich)
  plot(r_burn)
  plot(r_precip)
  plot(r_dcv)
}

get_values <- function(rst, pts) {
  #' extract values from `rst` after projecting `pts` to `crs(raster)`
  extract(rst, select(pts, x, y) %>%
            vect(geom = c('x', 'y')) %>%
            set.crs('EPSG:4326') %>%
            project(crs(rst))) %>%
    pull(2)
}

d <- as.data.frame(r_mu, xy = TRUE) %>%
  as_tibble() %>%
  mutate(.,
         s2_hat = get_values(r_s2, .),
         hfi = get_values(r_hfi, .), # ranges [0, 1]
         richness = get_values(r_rich, .),
         burned = get_values(r_burn, .),
         precip_m_year = get_values(r_precip, .) * 365,
         dcv_risk = get_values(r_dcv, .)) %>%
  mutate(mu_hat = if_else(mu_hat < -0.1, -0.1, mu_hat),
         s2_hat = if_else(s2_hat < 0, 0, s2_hat),
         s2_hat = if_else(s2_hat > 0.05, 0.05, s2_hat),
         precip_m_year = if_else(precip_m_year > 5, 5, precip_m_year)) %>%
  pivot_longer(- c(x, y, s2_hat), names_to = 'variable',
               values_to = 'value') %>%
  mutate(lab = case_when(
    variable == 'mu_hat' ~ 'Estimated mean NDVI',
    variable == 'hfi' ~ 'Human footprint',
    variable == 'richness' ~ 'Species richness',
    variable == 'dcv_risk' ~ 'Risk of arboviral diseases',
    variable == 'burned' ~ 'Mean proportion buned',
    variable == 'precip_m_year' ~ 'Annual precipitation (m)') %>%
      factor(., levels = unique(.)))
d

fig_4 <-
  ggplot() +
  facet_wrap(~ lab, scales = 'free_x', strip.position = 'bottom', nrow = 2) +
  geom_hex(aes(value, s2_hat, fill = log10(after_stat(count))), d,
           color = 'black', bins = 30, linewidth = 0.1, na.rm = TRUE) +
  scale_fill_lapaz(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, labels = \(.x) 10^.x) +
  labs(x = NULL) +
  scale_y_continuous('DENVar') +
  theme(strip.background = element_blank(), strip.placement = 'outside',
        legend.position = 'top', legend.key.width = rel(2),
        strip.text = element_text(size = rel(1)))

ggsave('figures/fig-4.png', fig_4,
       width = 12, height = 7, units = 'in', dpi = 600, bg = 'white')
