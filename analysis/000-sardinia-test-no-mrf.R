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
  st_geometry() %>%
  st_make_valid() %>%
  st_intersection(sardinia_bbox)

# import ndvi data (downloaded from the NASA website)
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
  saveRDS(sardinia_ndvi, 'data/sardinia-test/sardinia-ndvi.rds')
}

if(FALSE) {
  # plot the first few NDVI rasters
  filter(sardinia_ndvi, date %in% unique(date)[1:21]) %>%
    ggplot() +
    facet_wrap(~ date, nrow = 3) +
    coord_equal() +
    geom_raster(aes(x, y, fill = ndvi)) +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
}

# add altitude (get_elev_points results in a json error)
elevs <- sardinia_ndvi %>%
  filter(date == first(date)) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 6)
plot(elevs)
plot(sardinia, add = TRUE, col = 'transparent')

sardinia_ndvi <- mutate(sardinia_ndvi, 
                        elev_m = extract(elevs, select(sardinia_ndvi, x, y)),
                        year = year(date),
                        doy = yday(date) / 366) %>%
  filter(elev_m > -100)
sardinia_ndvi

# fit a model with gaussian family ----
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
draw(m_gaus, rug = FALSE, fun = \(x) (x + coef(m_gaus)['(Intercept)']) * 2 - 1)
ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-terms-w-intercept.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
summary(m_gaus)

# fit a test model with beta family ----
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
draw(m_beta, rug = FALSE, fun = \(x) m_beta$family$linkinv(x + coef(m_beta)['(Intercept)']) * 2 - 1)
ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms-ndvi-scale.png',
       width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
summary(m_beta)

# fit a test model with betals family ----
#' fails with error (tested on lab linux):
#'`Error in h(simpleError(msg, call)) :`
#'`error in evaluating the argument 'x' in selecting a method for function 't':`
#'`long vectors (argument 5) are not supported in .Fortran`
#'`Timing stopped at: 3.938e+04 5.971e+04 9.906e+04`
if(file.exists('models/sardinia-test/betals-gam.rds')) {
  m_betals <- readRDS('models/sardinia-test/betals-gam.rds')
} else {
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
  saveRDS(m_betals, 'models/sardinia-test/betals-gam.rds')
}

draw(m_betals, rug = FALSE)
summary(m_betals)

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
    scale_fill_gradientn(expression('NDVI,'~nu), colors = ndvi_pal, limits = c(-1, 1)) +
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
