library('dplyr')   # for data wrangling
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('sf')      # for simple fleature objects
library('terra')   # for rasters
library('gratia')  # for plotting GAMs
source('analysis/figures/000-default-ggplot-theme.R')

ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# import the data for the model
d <- readRDS('data/hbam-var-ndvi-data.rds')
nbs <- readRDS('data/ecoregions/poly-nbs-global.rds')

if(length(list.files('models/global-test', 'hbam-var-ndvi-*')) > 1) {
  DATE <- '2025-04-XX'
  m <- readRDS(paste0('models/global-test/hbam-var-ndvi-', DATE, '.rds'))
} else {
  m <- bam(
    e_2 ~
      biome + # to avoid intercept shrinkage
      s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
      s(doy, biome, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
      s(year, biome, bs = 'fs', xt = list(bs = 'cr'), k = 10) +
      ti(doy, year, biome, bs = c('cc', 'cr', 're'), k = c(5, 5)) +
      s(elevation_m, bs = 'cr', k = 5) +
      s(mu_hat, bs = 'cr', k = 5),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    drop.unused.levels = TRUE, # FALSE slows down the model substantially
    discrete = TRUE,
    samfrac = 0.001, # find intial guesses with a subset of the data
    nthreads = future::availableCores(logical = FALSE) - 2,
    control = gam.control(trace = TRUE))
  
  saveRDS(m, paste0('models/global-test/hbam-var-ndvi-', DATE, '.rds'))
}

# plot each smooth separately
p_hbam_doy <- draw(m, rug = FALSE,
                   select = which(grepl('s(doy', smooths(m))))
ggsave(paste0('figures/hbam-var-ndvi-fe-doy-', DATE, '.png'),
       plot = p_hbam_doy, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_y <- draw(m, rug = FALSE,
                 select = which(grepl('s(year', smooths(m))))
ggsave(paste0('figures/hbam-var-ndvi-fe-year-', DATE, '.png'),
       plot = p_hbam_y, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_elev <- draw(m, rug = FALSE,
                    select = which(smooths(m) == 's(elev_m)'))
ggsave(paste0('figures/hbam-var-ndvi-fe-elev_m-', DATE, '.png'),
       plot = p_hbam_elev, width = 4, height = 6, units = 'in', dpi = 300,
       bg = 'white')

p_hbam <- draw(m, rug = FALSE, - which(smooths(m) == 's(poly_id)'))

ggsave(paste0('figures/hbam-var-ndvi-', DATE, '.png'), plot = p_hbam,
       width = 8, height = 12, units = 'in', dpi = 600, bg = 'white')

# hex plot of mean and variance
ggplot() +
  geom_hex(aes(mu_hat, fitted(m))) +
  labs(x = 'Estimated mean NDVI', y = 'Estimated variance in NDVI')

ggsave('figures/hbam-var-mean-ndvi-hex.png',
       width = 8, height = 12, units = 'in', dpi = 600, bg = 'white')
