# hex plot of mean and variance
d %>%
  mutate(s2_hat = fitted(m_s2)) %>%
  filter(s2_hat < 0.25) %>% #' *remove*
  ggplot() +
  geom_hex(aes(mu_hat, s2_hat, fill = log10(after_stat(count))),
           color = 'black', bins = 100, linewidth = 0.1, na.rm = TRUE) +
  geom_hline(yintercept = 0, color = 'grey') +
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
  vect() %>% # convert to spatVect
  project(crs(dem))
m_mu <- readRDS('models/global-models/hbam-mean-ndvi-sos-5-no-res-2025-04-21-THINNED-50.rds')

preds <- dem %>%
  crop(ecoregions) %>%
  mask(ecoregions) %>%
  as.data.frame(dem, xy = TRUE) %>%
  rename(elevation_m = 3) %>%
  mutate(doy = 0, year = 0) %>%
  mutate(mu_hat = predict.bam(m_mu, newdata = ., discrete = FALSE,
                              terms = c('(Intercept)', 's(y,x)', 's(elevation_m)'))) %>%
  mutate(s2_hat = predict(m_s2, newdata = ., discrete = FALSE,
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
