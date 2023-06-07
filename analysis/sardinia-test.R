library('sf')
library('raster')
library('elevatr')
library('dplyr')
library('mgcv')
library('ggplot2')
source('functions/betals.r') # custom beta location-scale family
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
theme_set(theme_bw())

# shapefile of italy
italy <- filter(spData::world, name_long == 'Italy')

# import ndvi data (downloaded from the NASA website)
if(file.exists('data/sardinia-ndvi.rds')) {
  sardinia <- readRDS('data/sardinia-ndvi.rds')
} else {
  d <- readRDS('data/nasa-rasters/all-nasa-rasters.rds')
  sardinia <- filter(d,
                     long > 8, long < 10,
                     lat > 38, lat < 41.3) %>%
    filter(ndvi != 1) # crude way of dropping the water pixels
  rm(d)
  gc() # rm() doesn't always clear the RAM fully
  saveRDS(sardinia, 'data/sardinia-ndvi.rds')
}

# plot the first NDVI raster
filter(sardinia, date == min(date)) %>%
  ggplot() +
  coord_equal() +
  geom_raster(aes(long, lat, fill = ndvi)) +
  scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))

# some odd values < 0.5 in the south...
if(FALSE) {
  for(y in unique(sardinia$year)) {
    if(y > min(sardinia$year)) readline('Press <Enter> for the next plot.')
    (filter(sardinia, year == y) %>%
        ggplot() +
        facet_wrap(~ date, nrow = 3) +
        coord_equal() +
        geom_raster(aes(long, lat, fill = ndvi)) +
        ggtitle(y) +
        scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))) %>%
      print()
  }
}
# add altitude (get_elev_points results in a json error)
elevs <- sardinia %>%
  filter(date == min(date)) %>%
  select(long, lat) %>%
  SpatialPoints(proj4string = crs('EPSG:4326')) %>%
  raster() %>%
  get_elev_raster(z = 2)
plot(elevs)
plot(italy, add = TRUE, col = 'transparent')

sardinia$elevation <- extract(elevs, select(sardinia, long, lat))

# fit a test model with betals
m <- gam(list(ndvi_scaled ~ s(long, lat, bs = 'ds', k = 50) +
                s(elevation, bs = 'cr', k = 4) +
                s(year, bs = 'cr', k = 10) +
                s(doy, bs = 'cc', k = 10),
              ~ s(long, lat, bs = 'ds', k = 20) +
                s(elevation, bs = 'cr', k = 4) +
                s(year, bs = 'cr', k = 5) +
                s(doy, bs = 'cc', k = 5)),
         family = betals(),
         knots = list(doy = c(0, 1)),
         data = sardinia,
         method = 'REML',
         control = gam.control(nthreads = 11, trace = TRUE))

png('figures/sardinia-test/sardinia-ndvi-betals-terms.png', width = 8,
    height = 10, units = 'in', res = 300)
plot(m, pages = 1, scheme = 3, too.far = 0.05)
dev.off()

png('figures/sardinia-test/sardinia-ndvi-betals-gam-check.png', width = 8,
    height = 10, units = 'in', res = 300)
layout(matrix(1:4, ncol = 2))
gam.check(m, type = 'pearson')
dev.off()

# plots of predictions ----
newd <- filter(sardinia, date == min(date)) # predictions for first day

preds <-
  bind_cols(newd,
            as.data.frame(predict(m, newdata = newd, type = 'response')) %>%
              rename(mu = V1, # mean parameers
                     phi = V2) %>% # scale parameter
              mutate(sigma2 = phi * (1 - mu) * mu, # calculate variance
                     mu = mu * 2 - 1, # rescale to [-1, 1]
                     sigma2 = sigma2 * 4)) # scale variance appropriately

cowplot::plot_grid(
  # plot estimated mean
  ggplot() +
    coord_equal() +
    geom_raster(aes(long, lat, fill = mu), preds) +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1)),
  # plot estimated variance
  ggplot() +
    coord_equal() +
    geom_raster(aes(long, lat, fill = sigma2), preds) +
    scale_fill_viridis_c(limits = c(0, 1), values = seq(0, 1, by = 0.1)^2) # full range of sigma2
)

ggsave('figures/sardinia-test/sardinia-ndvi-betals-predictions.png',
       height = 10, width = 16)
