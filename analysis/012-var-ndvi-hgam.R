library('dplyr')   # for data wrangling
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('sf')      # for simple fleature objects
library('terra')   # for rasters
library('gratia')  # for plotting GAMs
library('ggplot2') # for fancy plots
source('analysis/figures/000-default-ggplot-theme.R')

# import the data for the model
d <- readRDS('data/hbam-var-ndvi-data-mod-5-no-res-2025-04-21-THINNED-50.rds')

hist(d$mu_hat)
quantile(d$mu_hat, 0.001)
# d <- d %>%
#   mutate(mu_hat = if_else(mu_hat < quantile(mu_hat, 0.001),
#                           quantile(mu_hat, 0.001), mu_hat),
#          mu_hat = if_else(mu_hat > 1, 1, mu_hat))
# nbs <- readRDS('data/ecoregions/poly-nbs-global.rds')

if(length(list.files('models/global-test', 'hbam-var-ndvi-*')) > 1) {
  DATE <- '2025-04-XX'
  m_s2 <- readRDS(paste0('models/global-models/hbam-var-ndvi-sos-mod-5-no-res-2025-04-21-THINNED-50.rds'))
} else {
  m_s2 <- bam(
    e_2 ~
      s(doy, bs = 'cc', k = 10) +
      s(year, bs = 'cr', k = 9) +
      s(y, x, bs = 'sos', k = 1e3) + # sos uses lat, long
      ti(doy, y, x, bs = c('cc', 'sos'), d = c(1, 2), k = c(10, 250)) +
      ti(year, y, x, bs = c('cr', 'sos'), d = c(1, 2), k = c(9, 250)) +
      ti(doy, year, y, x, bs = c('cc', 'cr', 'sos'), d = c(1, 1, 2),
         k = c(5, 5, 100)) +
      s(elevation_m, bs = 'cr', k = 5) +
      s(mu_hat, k = 5),
    family = gaussian(),
    data = d,
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    discrete = TRUE,
    samfrac = 0.01, # find initial guesses with a subset of the data
    nthreads = future::availableCores(logical = FALSE) - 2,
    control = gam.control(trace = TRUE))
  
  saveRDS(m_s2, paste0('models/global-models/hbam-var-ndvi-sos-', DATE, '.rds'))
}

100 * (1 - m_s2$deviance / m_s2$null.deviance)

png(paste0('figures/global-models/var-ndvi-terms-', DATE, '.png'),
    width = 12, height = 9, units = 'in', bg = 'white', res = 300)
plot(m_s2, pages = 1, too.far = 0.01, scheme = c(1, 1, 0, 5, 5, 5, 1, 1),
     phi = 90, theta = 0, n2 = 100, scale = 0)
dev.off()

# plot each smooth separately
p_hbam_doy <- draw(m_s2, rug = FALSE,
                   select = which(grepl('s(doy', smooths(m_s2))))
ggsave(paste0('figures/hbam-var-ndvi-fe-doy-', DATE, '.png'),
       plot = p_hbam_doy, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_y <- draw(m_s2, rug = FALSE,
                 select = which(grepl('s(year', smooths(m_s2))))
ggsave(paste0('figures/hbam-var-ndvi-fe-year-', DATE, '.png'),
       plot = p_hbam_y, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_elev <- draw(m_s2, rug = FALSE,
                    select = which(smooths(m_s2) == 's(elev_m)'))
ggsave(paste0('figures/hbam-var-ndvi-fe-elev_m-', DATE, '.png'),
       plot = p_hbam_elev, width = 4, height = 6, units = 'in', dpi = 300,
       bg = 'white')

p_hbam <- draw(m_s2, rug = FALSE, - which(smooths(m_s2) == 's(poly_id)'))

ggsave(paste0('figures/hbam-var-ndvi-', DATE, '.png'), plot = p_hbam,
       width = 8, height = 12, units = 'in', dpi = 600, bg = 'white')

# hex plot of mean and variance
d %>%
  mutate(s2_hat = fitted(m_s2)) %>%
  filter(s2_hat < 0.25) %>% #' *remove*
  ggplot() +
  geom_hex(aes(mu_hat, s2_hat, fill = log10(after_stat(count))),
           color = 'black', bins = 100, linewidth = 0.1, na.rm = TRUE) +
  labs(x = 'Estimated mean NDVI', y = 'DENVar') +
  scale_fill_iridescent(
    name = expression(paste(bold('Count (log'), bold(''['10']),
                            bold(' scale)'))), range = c(0, 1),
    reverse = FALSE, breaks = seq(0, 6, by = 1),
    labels = \(.x) 10^.x) +
  theme(legend.position = 'top', legend.key.width = unit(0.5, 'in'))

ggsave(paste0('figures/global-models/hbam-var-mean-ndvi-hex-', DATE, '.png'),
       width = 8, height = 6, units = 'in', dpi = 600, bg = 'white')

# make rasters of space
dem <- rast('data/elev-raster.tif')
ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  # drop polygons not in data
  filter(poly_id %in% names(readRDS('data/ecoregions/poly-nbs-global.rds'))) %>%
  vect() %>% # convert to spatVect
  project(crs(dem))

preds <- dem %>%
  crop(ecoregions) %>%
  mask(ecoregions) %>%
  as.data.frame(dem, xy = TRUE) %>%
  rename(elevation_m = 3) %>%
  mutate(doy = 0, year = 0) %>%
  mutate(mu_hat = predict.bam(m_mu, newdata = ., discrete = TRUE,
                          terms = c('(Intercept)', 's(y,x)', 's(elevation_m)'))) %>%
  mutate(s2_hat = predict(m_s2, newdata = ., discrete = TRUE,
                          terms = c('(Intercept)', 's(y,x)', 's(elevation_m)')))
head(preds)

r_mu <- preds %>%
  select(x, y, mu_hat) %>%
  rast(crs = crs(dem))

r_s2 <- preds %>%
  select(x, y, s2_hat) %>%
  rast(crs = crs(dem))

plot(r_mu)
plot(r_s2)
plot(values(r_s2) ~ values(r_mu))

writeRaster(r_mu, paste0('data/output/mean-ndvi-raster-', DATE, '.tif'), overwrite = TRUE)
writeRaster(r_s2, paste0('data/output/var-ndvi-raster-', DATE, '.tif'), overwrite = TRUE)

rast('H:/GitHub/ndvi-stochasticity/data/output/mean-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif') %>%
  as.data.frame(xy = TRUE) %>%
  mutate(mu_hat = case_when(
    mu_hat < quantile(mu_hat, 0.0001) ~ quantile(mu_hat, 0.0001),
    mu_hat > quantile(mu_hat, 0.9999) ~ quantile(mu_hat, 0.9999),
    .default = mu_hat)) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = mu_hat)) +
  geom_sf(data = st_as_sf(ecoregions), color = 'grey30', linewidth = 0.1,
          fill = 'transparent', inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c('Mean NDVI') +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = unit(0.5, 'in'))
ggsave('H:/GitHub/ndvi-stochasticity/figures/global-models/output/mean-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.png',
       width = 8, height = 4, units = 'in', bg = 'white', dpi = 300)

rast('H:/GitHub/ndvi-stochasticity/data/output/var-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.tif') %>%
  as.data.frame(xy = TRUE) %>%
  mutate(s2_hat = case_when(
    s2_hat < quantile(s2_hat, 0.0001) ~ quantile(s2_hat, 0.0001),
    s2_hat > 0.04 ~ 0.04,
    .default = s2_hat)) %>%
  ggplot() +
  geom_raster(aes(x, y, fill = s2_hat)) +
  geom_sf(data = st_as_sf(ecoregions), color = 'grey30', linewidth = 0.1,
          fill = 'transparent', inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_viridis_c('DENVar, capped at 0.04', option = 'A') +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = unit(0.5, 'in'))
ggsave('H:/GitHub/ndvi-stochasticity/figures/global-models/output/var-ndvi-raster-mod-5-no-res-2025-04-21-THINNED-50.png',
       width = 8, height = 4, units = 'in', bg = 'white', dpi = 300)
