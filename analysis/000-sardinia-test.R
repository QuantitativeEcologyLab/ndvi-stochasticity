# all tests ran on my personal laptop:
# 64.0 GB RAM, 13th Gen Intel Core i7-1370P processor (1.90 GHz)
library('sf')
library('terra')
library('elevatr')
library('dplyr')
library('lubridate')
library('purrr')
library('furrr')
library('mgcv')
library('ggplot2')
library('gratia')
source('functions/betals.r') # custom beta location-scale family
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('analysis/figures/default-ggplot-theme.R')
theme_set(theme_bw())

# the bounding box for sardinia (large enough to include all coast)
sardinia_bbox <- tibble(x = c(7.6, 10.4), y = c(38.8, 41.3)) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_bbox() %>%
  st_as_sfc() %>%
  st_as_sf() %>%
  st_set_crs('EPSG:4326')

# shapefile of the world's terrestrial ecoregions
sardinia <- read_sf('data/world-ecosystems/data/commondata/data0/tnc_terr_ecoregions.shp') %>%
  st_make_valid() %>%
  st_intersection(sardinia_bbox) %>%
  st_cast('POLYGON') %>%
  st_as_sf() %>%
  mutate(poly_id = paste('poly', 1:n()),
         area_km2 = as.numeric(st_area(.)) / 1e6)
hist(sardinia$area_km2)

# import ndvi data (downloaded from the NASA website) ----
if(file.exists('data/sardinia-test/sardinia-ndvi.rds')) {
  sardinia_ndvi <- readRDS('data/sardinia-test/sardinia-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix') stop('NOAA rasters are on the lab Linux.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = 30)
  sardinia_ndvi <-
    list.files(path = '/home/shared/NOAA_Files/',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn)
      r %>%
        crop(sardinia_bbox) %>% # to avoid including the rest of italy
        mask(sardinia_bbox) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(ndvi_scaled = ndvi_to_01(ndvi))
  plan(sequential)
  
  # plot the first few NDVI rasters
  filter(sardinia_ndvi, date %in% unique(date)[1:21]) %>%
    ggplot() +
    facet_wrap(~ date, nrow = 3) +
    coord_equal() +
    geom_raster(aes(x, y, fill = ndvi)) +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  # add altitude (get_elev_points results in a json error)
  elevs <- sardinia_ndvi %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 6)
  plot(elevs)
  plot(st_geometry(sardinia), add = TRUE, col = 'transparent')
  
  sardinia_ndvi <- mutate(sardinia_ndvi, 
                          elev_m = extract(elevs, select(sardinia_ndvi, x, y)),
                          year = year(date),
                          doy = yday(date) / 366) %>%
    filter(elev_m > -100)
  sardinia_ndvi
  
  saveRDS(sardinia_ndvi, 'data/sardinia-test/sardinia-ndvi.rds')
}

# fit a spatially explicit test model with a gaussian family ----
# on personal laptop: "initial" is 26 s, run is ~2 s
if(file.exists('models/sardinia-test/gaussian-gam.rds')) {
  m_gaus <- readRDS('models/sardinia-test/gaussian-gam.rds')
} else {
  m_gaus <- bam(ndvi_scaled ~ s(x, y, bs = 'ds', k = 200) +
                  s(elev_m, bs = 'ad', k = 10) + # adaptive spline to correct for shorelines better
                  s(year, bs = 'cr', k = 10) +
                  s(doy, bs = 'cc', k = 10),
                family = gaussian(),
                knots = list(doy = c(0, 1)),
                data = sardinia_ndvi,
                method = 'fREML',
                discrete = TRUE,
                control = gam.control(nthreads = 1, trace = TRUE))
  saveRDS(m_gaus, 'models/sardinia-test/gaussian-gam.rds')
}

draw(m_gaus, rug = FALSE)
ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-terms.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
draw(m_gaus, rug = FALSE, fun = \(x) ndvi_to_11(x + coef(m_gaus)['(Intercept)']))
ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-terms-ndvi-scale.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
summary(m_gaus)

# fit a spatially explicit test model with a beta family ----
# on personal laptop: "initial" is 30 s, run is ~ 660 s
if(file.exists('models/sardinia-test/beta-gam.rds')) {
  m_beta <- readRDS('models/sardinia-test/beta-gam.rds')
} else {
  m_beta <- bam(ndvi_scaled ~ s(x, y, bs = 'ds', k = 200) +
                  s(elev_m, bs = 'ad', k = 10) +
                  s(year, bs = 'cr', k = 10) +
                  s(doy, bs = 'cc', k = 10),
                family = betar(),
                knots = list(doy = c(0, 1)),
                data = sardinia_ndvi,
                method = 'fREML',
                discrete = TRUE,
                control = gam.control(trace = TRUE))
  saveRDS(m_beta, 'models/sardinia-test/beta-gam.rds')
}

draw(m_beta, rug = FALSE)
ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
draw(m_beta, rug = FALSE,
     fun = \(x) m_beta$family$linkinv(x + coef(m_beta)['(Intercept)']) %>%
       ndvi_to_11())
ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms-ndvi-scale.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
summary(m_beta)

# fit a spatially explicit test model with a betals family ----
#' fails with error (*tested on lab linux machine*):
#'`Error in h(simpleError(msg, call)) :`
#'`error in evaluating the argument 'x' in selecting a method for function 't':`
#'`long vectors (argument 5) are not supported in .Fortran`
#'`Timing stopped at: 3.938e+04 5.971e+04 9.906e+04`
if(FALSE) {
  system.time(
    m_betals <- gam(list(ndvi_scaled ~ s(x, y, bs = 'ds', k = 200) +
                           s(elev_m, bs = 'ad', k = 10) +
                           s(year, bs = 'cr', k = 10) +
                           s(doy, bs = 'cc', k = 10),
                         ~ s(x, y, bs = 'ds', k = 200) +
                           s(elev_m, bs = 'ad', k = 10) +
                           s(year, bs = 'cr', k = 10) +
                           s(doy, bs = 'cc', k = 10)),
                    family = betals(),
                    knots = list(doy = c(0, 1)),
                    data = sardinia_ndvi,
                    method = 'REML',
                    control = gam.control(nthreads = 10, trace = TRUE)))
}

#' fitting gaussian models rather than beta models because they are
#' substantially faster with no clear loss to predictive accuracy
p_fits <-
  tibble(beta = ndvi_to_11(fitted(m_beta)),
         gaus = ndvi_to_11(fitted(m_gaus))) %>%
  ggplot() +
  geom_point(aes(x, x), data = tibble(x = 0:1), color = 'transparent') +
  geom_point(aes(beta, gaus),
             alpha = 0.01) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Beta model (with constant scale parameter)',
       y = 'NDVI predicted by Gaussian model (with constant variance)')
ggsave('sardinia-ndvi-model-agreement-gaussian-beta.png', plot = p_fits,
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

# plots of predicted means and squared residuals for day 1 ----
preds_1 <- sardinia_ndvi %>%
  mutate(mu_gaus = ndvi_to_11(fitted(m_gaus)),
         mu_beta = ndvi_to_11(fitted(m_beta))) %>%
  filter(date == min(date))

ggplot() +
  geom_point(aes(mu_beta, mu_gaus), data = preds_1, alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Beta model (with constant scale parameter)',
       y = 'NDVI predicted by Gaussian model (with constant variance)')
ggsave(paste0('sardinia-ndvi-model-agreement-gaussian-beta-',
              unique(preds_1$date), '.png'),
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

preds_1_long <- preds_1 %>%
  tidyr::pivot_longer(cols = c(mu_gaus, mu_beta), names_to = 'model',
                      names_prefix = 'mu_', values_to = 'mu') %>%
  mutate(model = case_when(model == 'beta' ~ 'Beta',
                           model == 'gaus' ~ 'Gaussian'),
         e2 = (ndvi - mu)^2)

cowplot::plot_grid(
  # plot estimated mean
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = mu), preds_1_long) +
    scale_fill_viridis_c(expression('NDVI,'~nu), option = 'A') +
    labs(title = paste('Predictions for', sardinia_ndvi$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_long) +
    scale_fill_viridis_c(expression((nu-E(nu))^2), limits = c(0, NA)) +
    labs(title = paste('Squared residuals', sardinia_ndvi$date[1]), x = NULL, y = NULL),
  labels = 'AUTO', ncol = 1)

ggsave('figures/sardinia-test/sardinia-ndvi-model-agreement-gaussian-beta-predictions.png',
       height = 12, width = 8, bg = 'white')

# fit a gaussian model with markov random fields ----
ggplot(sardinia) + geom_sf()

# could use matrices of coordinates instead of list of neighbors, but not
# using this because some polygons wrap around the globe, so the nbs need
# to be fixed manually, and doing so with matrices would not work well 
st_coordinates(sardinia) %>%
  as.data.frame() %>%
  as_tibble() %>%
  select(! L1) %>%
  tidyr::nest(poly = ! L2) %>%
  mutate(poly = map(poly, as.matrix)) %>%
  pull(poly) %>%
  `names<-`(1:6)

# using polygons and neighbors (~ 10 seconds)
sardinia_ndvi <-
  sardinia_ndvi %>%
  tidyr::nest(data = ! c(x, y)) %>%
  mutate(poly_i = map2_int(x, y, function(long, lat) {
    .p <- st_as_sf(tibble(long, lat), coords = c('long', 'lat')) %>%
      st_set_crs('EPSG:4326') %>%
      st_intersects(sardinia, sparse = TRUE)
    
    # some raster cells are assigned to 0 or > 2 polygons (rarely)
    if(length(.p) == 0) {
      return(NA_integer_)
    } else if(length(.p) == 1) {
      return(as.integer(.p))
    } else {
      warning('Selected the first of ', length(.p), ' polygons.')
      return(as.integer(.p)[1])
    }
  })) %>%
  filter(! is.na(poly_i)) %>%
  mutate(eco_name = factor(sardinia$ECO_NAME[poly_i]),
         poly_id = sardinia$poly_id[poly_i],
         poly_id = factor(poly_id, levels = unique(sardinia$poly_id))) %>%
  select(! poly_i) %>%
  tidyr::unnest(data)

# find neighbors for each polygon
nbs <-
  spdep::poly2nb(pl = sardinia,
                 row.names = sardinia$poly_id, # so neighbors match ID
                 queen = TRUE) # 1 point in common is sufficient

# no neighbors because all polygons are islands
nbs[1:length(nbs)] # ith element is the integer index of the neighbors
length(nbs) == nrow(sardinia)

# add names corresponding to each polygon
all(attr(nbs, 'region.id') == sardinia$poly_id)
names(nbs) # no names, currently
names(nbs) <- attr(nbs, 'region.id') # names need to match names
names(nbs) # polygons' names (same as sardinia_ndvi$poly_id)

if(file.exists('models/sardinia-test/gaussian-mrf-gam.rds')) {
  m_mrf <- readRDS('models/sardinia-test/gaussian-mrf-gam.rds')
} else {
  # mrf only
  m_mrf_0 <- bam(ndvi_scaled ~ s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
                 family = gaussian(),
                 data = sardinia_ndvi,
                 method = 'fREML',
                 discrete = TRUE,
                 control = gam.control(trace = TRUE))
  
  # mrf, elevation, and time
  m_mrf_1 <- bam(ndvi_scaled ~
                   s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
                   s(elev_m, bs = 'ad', k = 10) + # adaptive spline to correct for shorelines better
                   s(year, bs = 'cr', k = 10) +
                   s(doy, bs = 'cc', k = 10),
                 family = gaussian(),
                 knots = list(doy = c(0, 1)),
                 data = sardinia_ndvi,
                 method = 'fREML',
                 discrete = TRUE,
                 control = gam.control(trace = TRUE))
  
  # adding a complex smooth of explicit space
  m_mrf_2 <- bam(ndvi_scaled ~
                   s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
                   s(x, y, bs = 'ds', k = 200) +
                   s(elev_m, bs = 'ad', k = 10) +
                   s(year, bs = 'cr', k = 10) +
                   s(doy, bs = 'cc', k = 10),
                 family = gaussian(),
                 knots = list(doy = c(0, 1)),
                 data = sardinia_ndvi,
                 method = 'fREML',
                 discrete = TRUE,
                 control = gam.control(trace = TRUE))
  
  round(summary(m_mrf_0)$dev.expl * 100, 2)
  round(summary(m_mrf_1)$dev.expl * 100, 2) # increases substantially
  round(summary(m_mrf_2)$dev.expl * 100, 2) # does not increase much
  
  saveRDS(m_mrf_1, 'models/sardinia-test/gaussian-mrf-gam.rds')
  
  if(FALSE) {
    # even a simple beta model with only mrfs is much slower
    # fits in 19.7 minutes (compared to < 5 seconds for the gaussian model)
    m_mrf_b <- bam(ndvi_scaled ~ s(poly_id, bs = 'mrf', xt = list(nb = nbs)),
                   family = betar(),
                   data = sardinia_ndvi,
                   method = 'fREML',
                   discrete = TRUE,
                   control = gam.control(trace = TRUE))
  }
  
  m_mrf <- m_mrf_1
  rm(m_mrf_0, m_mrf_1, m_mrf_2, m_mrf_b)
}

# can't draw mrf smooths fit using polygons
draw(m_mrf, rug = FALSE, select = -1, scales = 'fixed')
summary(m_mrf)

mrf_preds <- sardinia %>%
  mutate(year = 0, doy = 0, elev_m = 0) %>%
  mutate(mu = predict(m_mrf, ., type = 'response', terms = 's(poly_id)'),
         ndvi = ndvi_to_11(mu + coef(m_mrf)['(Intercept)']))

# recreate the plot from draw()
mrf_plots <- list(
  ggplot(mrf_preds) +
    geom_sf(aes(fill = mu)) +
    scale_fill_distiller('Partial\neffect', type = 'div', palette = 5,
                         limits = c(-1, 1) * max(abs(mrf_preds$mu))),
  draw(m_mrf, rug = FALSE, select = 2),
  draw(m_mrf, rug = FALSE, select = 3),
  draw(m_mrf, rug = FALSE, select = 4))
cowplot::plot_grid(plotlist = mrf_plots)
ggsave('figures/sardinia-test/sardinia-ndvi-mrf-terms.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')

mrf_plots_11 <- list(
  ggplot(mrf_preds) +
    geom_sf(aes(fill = ndvi)) +
    scale_fill_distiller('Partial\neffect', type = 'div', palette = 5,
                         limits = c(-1, 1) * max(abs(mrf_preds$ndvi))),
  draw(m_mrf, rug = FALSE, select = 2,
       fun = \(x) ndvi_to_11(x + coef(m_mrf)['(Intercept)'])),
  draw(m_mrf, rug = FALSE, select = 3,
       fun = \(x) ndvi_to_11(x + coef(m_mrf)['(Intercept)'])),
  draw(m_mrf, rug = FALSE, select = 4,
       fun = \(x) ndvi_to_11(x + coef(m_mrf)['(Intercept)'])))
cowplot::plot_grid(plotlist = mrf_plots_11)
ggsave('figures/sardinia-test/sardinia-ndvi-mrf-terms-ndvi-scale.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')

# relatively close deviance explained for the three models
# (considering that the deviances are different between gaussian and beta)
round(summary(m_gaus)$dev.expl * 100, 2)
round(summary(m_beta)$dev.expl * 100, 2)
round(summary(m_mrf)$dev.expl * 100, 2)

# check agreement with other models
p_fits_mrf <-
  tibble(ds = predict(m_gaus, newdata = sardinia_ndvi, type = 'response') %>%
           ndvi_to_11(),
         mrf = predict(m_mrf, newdata = sardinia_ndvi, type = 'response') %>%
           ndvi_to_11()) %>%
  ggplot() +
  geom_point(aes(x, x), data = tibble(x = 0:1), color = 'transparent') +
  geom_point(aes(ds, mrf), alpha = 0.01) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Gaussian model (with Duchon spatial smooth)',
       y = 'NDVI predicted by Gaussian model (with MRF smooth)')
ggsave('sardinia-ndvi-model-agreement-duchon-mrf.png', plot = p_fits_mrf,
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

#' to fix the discordance above, we could add the average residual for a
#' given pixel to the mean and subtract it from the residuals before
#' calculating the variance
#' 
#' could also correct for NDVI bias in islands using log(area_km2)

# plots of predicted means and squared residuals for day 1 ----
preds_1_mrf <- sardinia_ndvi %>%
  filter(date == min(date)) %>%
  mutate(.,
         ds = ndvi_to_11(predict(m_gaus, newdata = .,
                                 type = 'response', discrete = TRUE)),
         mrf = ndvi_to_11(predict(m_mrf, newdata = .,
                                  type = 'response', discrete = TRUE)))

ggplot() +
  geom_point(aes(ds, mrf), data = preds_1_mrf, alpha = 0.1) +
  geom_abline(intercept = 0, slope = 1, color = 'red') +
  labs(x = 'NDVI predicted by Gaussian model (with Duchon spatial smooth)',
       y = 'NDVI predicted by Gaussian model (with MRF smooth)')
ggsave(paste0('sardinia-ndvi-model-agreement-duchon-mrf-',
              unique(preds_1_mrf$date), '.png'),
       path = 'figures/sardinia-test', width = 9, height = 6, units = 'in',
       dpi = 300, bg = 'white')

preds_1_mrf_long <- preds_1_mrf %>%
  tidyr::pivot_longer(cols = c(ds, mrf), names_to = 'model', values_to = 'mu') %>%
  mutate(model = case_when(model == 'ds' ~ 'Duchon smooth',
                           model == 'mrf' ~ 'Markov Random Fields'),
         e2 = (ndvi - mu)^2)

cowplot::plot_grid(
  # plot estimated mean
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = mu), preds_1_mrf_long) +
    scale_fill_viridis_c(expression('NDVI,'~nu), option = 'A') +
    labs(title = paste('Predictions for', sardinia_ndvi$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_mrf_long) +
    scale_fill_viridis_c(expression((nu-E(nu))^2), limits = c(0, NA)) +
    labs(title = paste('Squared residuals', sardinia_ndvi$date[1]), x = NULL, y = NULL),
  labels = 'AUTO', ncol = 1)

ggsave(paste0('figures/sardinia-test/sardinia-ndvi-model-agreement-duchon-mrf-',
              unique(preds_1_mrf$date), 'predictions.png'),
       height = 12, width = 8, bg = 'white')

# testing complexity of ti(doy, space) ----
#' `ti()` of doy and space doesn't add much
# on personal laptop: "initial" is 50 s, run is 5 s
if(file.exists('models/sardinia-test/gaus-gam-ti-ds-gam.rds')) {
  m_gaus_ti_ds <- readRDS('models/sardinia-test/gaus-gam-ti-ds-gam.rds')
} else {
  m_gaus_ti_ds <- bam(ndvi_scaled ~
                      s(x, y, bs = 'ds', k = 200) +
                      s(elev_m, bs = 'ad', k = 10) + # adaptive spline to correct for shorelines better
                      s(year, bs = 'cr', k = 10) +
                      s(doy, bs = 'cc', k = 10) +
                      ti(doy, year, bs = c('cc', 'cr'), k = c(5, 5)),
                    family = gaussian(),
                    knots = list(doy = c(0, 1)),
                    data = sardinia_ndvi,
                    method = 'fREML',
                    discrete = TRUE,
                    control = gam.control(nthreads = 1, trace = TRUE))
  saveRDS(m_gaus_ti_ds, 'models/sardinia-test/gaus-gam-ti-ds-gam.rds')
}

draw(m_gaus_ti_ds, rug = FALSE)
ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-ti-terms.png',
       width = 13.5, height = 6, units = 'in', dpi = 300, bg = 'white')

summary(m_gaus_ti_ds)
