library('sf')      # for simple features
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for figures
library('khroma')  # for colorblind-friendly color scale
library('cowplot') # for ggplot plots in grids
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R')

r_s2 <- rast('output/long-term-preds.tif')[[2]]
ecoregions <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  vect() %>% # convert to spatVect
  project(crs(r_s2))
r_s2 <- crop(r_s2, ecoregions, mask = TRUE)
r_mu <- rast('output/long-term-preds.tif')[[1]]
r_dhi <- rast('data/other-rasters/dhi-data/dhi_ndvi_2015.tif') %>%
  project(r_s2)
#' **USING THE KEYS ONE FOR TESTS BECAUSE IT'S FASTER**
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

#' *data that might be worth adding:*
#' - temperature: `https://cds.climate.copernicus.eu/datasets/derived-near-surface-meteorological-variables?tab=download` (need to make a raster)
#' - seasonal temperature range (need to make a raster)
#' - max body size (can't find a raster)
#' - invasive species (can't find a raster)
#' - disease occurrence? (can't find a raster)

# some have a difference CRS, but reprojecting the HFI layer takes too long
tibble(raster = c('r_s2', 'r_mu', 'r_dhi', 'r_hfi', 'r_rich', 'r_burn'),
       crs = sapply(raster, \(x) crs(get(x), proj = TRUE)),
       same_crs = crs == crs[1])

# check rasters
if(FALSE) {
  plot(r_s2)
  plot(r_mu)
  plot(r_dhi)
  plot(r_hfi)
  plot(r_rich)
  plot(r_burn)
  plot(r_precip)
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
  mutate(.,
         s2_hat = get_values(r_s2, .),
         hfi = get_values(r_hfi, .) / 1e3, # scale back to [0, 50]
         richness = get_values(r_rich, .),
         burned = get_values(r_burn, .),
         precip_m_year = get_values(r_precip, .) * 365) %>%
  # DHI rasters need to be extracted separately
  bind_cols(.,
            extract(r_dhi, select(., x, y) %>%
                      vect(geom = c('x', 'y')) %>%
                      set.crs('EPSG:4326') %>%
                      project(crs(r_dhi))) %>%
              select(! ID) %>%
              rename(dhi_cumulative = dhi_ndvi_2015_1,
                     dhi_min = dhi_ndvi_2015_2,
                     dhi_seasonal = dhi_ndvi_2015_3)) %>%
  as_tibble()
d

fig_3 <-
  d %>%
  mutate(precip_m_year = if_else(precip_m_year > 5, 5, precip_m_year)) %>%
  tidyr::pivot_longer(- c(x, y, s2_hat), names_to = 'variable',
                      values_to = 'value') %>%
  mutate(lab = case_when(
    variable == 'mu_hat' ~ 'Estimated mean NDVI',
    variable == 'hfi' ~ 'Human footprint',
    variable == 'richness' ~ 'Species richness',
    variable == 'dhi_cumulative' ~ 'Cumulative NDVI',
    variable == 'dhi_min' ~ 'Minimum NDVI',
    variable == 'dhi_seasonal' ~ 'Seasonal CV of NDVI',
    variable == 'burned' ~ 'Mean proportion buned',
    variable == 'precip_m_year' ~ 'Annual precipitation (m)') %>%
      factor(., levels = unique(.))) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', strip.position = 'left', nrow = 2) +
  geom_hex(aes(s2_hat, value, fill = log10(after_stat(count))),
           color = 'black', bins = 30, linewidth = 0.1, na.rm = TRUE) +
  scale_fill_lapaz(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, labels = \(.x) 10^.x) +
  labs(y = NULL, x = 'DENVar') +
  theme(strip.background = element_blank(), strip.placement = 'outside',
        legend.position = 'top', legend.key.width = rel(2),
        strip.text = element_text(size = rel(1)))

ggsave('figures/fig-3.png', fig_3,
       width = 12, height = 5, units = 'in', dpi = 600, bg = 'white')
