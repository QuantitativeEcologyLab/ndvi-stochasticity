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
r_hfi <- rast('H:/GitHub/ndvi-stochasticity/data/other-rasters/hfp_2021_100m_v1-2_cog.tif')
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

# some rasters have a difference CRS, but reprojecting takes too long
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
         burned = get_values(r_burn, .)) %>%
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

if(FALSE) {
  d %>%
    filter(is.na(hfi)) %>%
    plot(y ~ x, .)
}

ggplot(d, aes(s2_hat)) +
  geom_histogram(fill = 'grey', color = 'black', binwidth = 0.005,
                 center = 0.0025) +
  labs(x = 'DENVar', y = 'Count')
ggsave('figures/denvar-hist-global-test.png',
       width = 5, height = 3, units = 'in', dpi = 600, bg = 'white')

d %>%
  filter(mu_hat > -0.25) %>% #' *remove: there were 4 points mean < -0.25*
  tidyr::pivot_longer(- c(x, y, s2_hat), names_to = 'variable',
                      values_to = 'value') %>%
  mutate(lab = case_when(variable == 'mu_hat' ~ 'Estimated mean NDVI',
                         variable == 'hfi' ~ 'Human footprint',
                         variable == 'richness' ~ 'Species richness',
                         variable == 'dhi_cumulative' ~ 'Cumulative NDVI',
                         variable == 'dhi_min' ~ 'Minimum NDVI',
                         variable == 'dhi_seasonal' ~ 'Seasonal CV of NDVI',
                         variable == 'burned' ~ 'Mean proportion buned') %>%
           factor(., levels = unique(.))) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free_y', strip.position = 'left', nrow = 2) +
  geom_hex(aes(s2_hat, value, fill = log10(after_stat(count))),
           color = 'black', bins = 50, linewidth = 0.1, na.rm = TRUE) +
  scale_fill_iridescent(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, labels = \(.x) 10^.x) +
  labs(y = NULL, x = 'DENVar') +
  theme(legend.position = 'inside', legend.position.inside = c(7/8, 0.25),
        strip.background = element_blank(),
        strip.placement = 'outside',# legend.key.width = unit(0.5, 'in'),
        strip.text = element_text(size = rel(1)))

ggsave('figures/hexplots-global-test.png',
       width = 12, height = 5, units = 'in', dpi = 600, bg = 'white')
