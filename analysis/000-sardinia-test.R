# all tests ran on my personal laptop:
# 64.0 GB RAM, 13th Gen Intel Core i7-1370P processor (1.90 GHz)
library('sf')        # for shapefiles
library('terra')     # for rasters
library('elevatr')   # for digital elevation models
library('dplyr')     # for data wrangling
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('mgcv')      # for GAMs
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids 
library('gratia')    # for fancy plots of GAMs
library('inlabru')   # for plotting mesh objects in ggplot
library('ggplot2')   # for fanct plots
library('fmesher')   # to generate triangle meshes in polygons
source('functions/betals.r') # custom beta location-scale family
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('analysis/figures/000-default-ggplot-theme.R')

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

sardinia_geom <- st_geometry(sardinia)

# import ndvi data (downloaded from the NASA website) ----
if(file.exists('data/sardinia-test/sardinia-ndvi.rds')) {
  d <- readRDS('data/sardinia-test/sardinia-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix') stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use nultiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(30, availableCores(logical = FALSE) - 2))
  d <-
    list.files(path = '/home/shared/NOAA_Files/',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI')
      r %>%
        crop(sardinia_geom) %>% # to avoid including the rest of italy
        mask(sardinia_geom) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    mutate(ndvi_scaled = ndvi_to_01(ndvi))
  plan(sequential)
  
  # plot the first few NDVI rasters
  filter(d, date %in% unique(date)[1:21]) %>%
    ggplot() +
    facet_wrap(~ date, nrow = 3) +
    coord_equal() +
    geom_raster(aes(x, y, fill = ndvi)) +
    scale_fill_gradientn(colors = ndvi_pal, limits = c(-1, 1))
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 6)
  plot(elevs)
  plot(sardinia_geom, add = TRUE, col = 'transparent')
  
  d <- mutate(d, 
                          elev_m = extract(elevs, select(d, x, y)),
                          year = year(date),
                          doy = yday(date)) %>%
    filter(elev_m > -100)
  
  # drop points outside sardinia
  # unique coordinates inside sardinia
  locs <- d %>%
    select(x, y) %>%
    group_by(x, y) %>%
    slice(1) %>%
    st_as_sf(coords = c('x', 'y')) %>%
    st_set_crs('EPSG:4326') %>%
    filter(., st_as_sf(., coords = c('x', 'y')) %>%
             st_set_crs('EPSG:4326') %>%
             st_intersects(sardinia_geom, sparse = TRUE) %>%
             map_lgl(\(x) length(x) > 0))
  nrow(locs) # all rasters use same coords
  
  plot(sardinia_geom)
  plot(locs, add = TRUE)
  
  d <- filter(d, paste(x, y) %in%
                            paste(st_coordinates(locs)[, 1],
                                  st_coordinates(locs)[, 2]))
  d
  
  plot(sardinia_geom)
  plot(locs, add = TRUE)
  
  saveRDS(d, 'data/sardinia-test/sardinia-ndvi.rds')
}
summary(d)

# unique coordinates inside sardinia
locs <- d %>%
  select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  filter(., st_as_sf(., coords = c('x', 'y')) %>%
           st_set_crs('EPSG:4326') %>%
           st_intersects(sardinia_geom, sparse = TRUE) %>%
           map_lgl(\(x) length(x) > 0))
nrow(locs) # all rasters use same coords

plot(sardinia_geom)
plot(locs, add = TRUE)

# fit a spatially explicit test model with a gaussian family ----
# on personal laptop: "initial" is 56 s, run is ~2 s
if(file.exists('models/sardinia-test/gaussian-gam.rds')) {
  m_gaus <- readRDS('models/sardinia-test/gaussian-gam.rds')
} else {
  m_gaus <- bam(ndvi_scaled ~ s(x, y, bs = 'ds', k = 200) +
                  s(elev_m, bs = 'ad', k = 10) + # adaptive spline to correct for shorelines better
                  s(year, bs = 'cr', k = 10) +
                  s(doy, bs = 'cc', k = 10),
                family = gaussian(),
                knots = list(doy = c(0, 1)),
                data = d,
                method = 'fREML',
                discrete = TRUE,
                control = gam.control(nthreads = 1, trace = TRUE))
  saveRDS(m_gaus, 'models/sardinia-test/gaussian-gam.rds')
  
  draw(m_gaus, rug = FALSE)
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  draw(m_gaus, rug = FALSE, fun = \(x) ndvi_to_11(x + coef(m_gaus)['(Intercept)']), dist = 0.03) &
    scale_fill_viridis_c('NDVI', limits = c(0, 0.4))
  ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-terms-ndvi-scale.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  summary(m_gaus)
}

# fit a spatially explicit test model with a beta family ----
# on personal laptop: "initial" is 35 s, run is ~ 25 minutes
if(file.exists('models/sardinia-test/beta-gam.rds')) {
  m_beta <- readRDS('models/sardinia-test/beta-gam.rds')
} else {
  m_beta <- bam(ndvi_scaled ~ s(x, y, bs = 'ds', k = 200) +
                  s(elev_m, bs = 'ad', k = 10) +
                  s(year, bs = 'cr', k = 10) +
                  s(doy, bs = 'cc', k = 10),
                family = betar(),
                knots = list(doy = c(0, 1)),
                data = d,
                method = 'fREML',
                discrete = TRUE,
                control = gam.control(trace = TRUE))
  saveRDS(m_beta, 'models/sardinia-test/beta-gam.rds')
  
  draw(m_beta, rug = FALSE)
  ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  draw(m_beta, rug = FALSE, fun = \(x) m_beta$family$linkinv(x + coef(m_beta)['(Intercept)']) %>%
         ndvi_to_11(), dist = 0.03) &
    scale_fill_viridis_c('NDVI', limits = c(0, 0.4))
  ggsave('figures/sardinia-test/sardinia-ndvi-beta-terms-ndvi-scale.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  summary(m_beta)
}

# fit a spatially explicit test model with a betals family ----
#' fails with error (*tested on lab linux machine*):
#'`Error in h(simpleError(msg, call)) :`
#'`error in evaluating the argument 'x' in selecting a method for function 't':`
#'`long vectors (argument 5) are not supported in .Fortran`
#'`Timing stopped at: 3.938e+04 5.971e+04 9.906e+04`
if(FALSE) {
  system.time(
    m_betals <- gam(list(ndvi_scaled ~ s(x, y, bs = 'ds', k = 20) +
                           s(elev_m, bs = 'ad', k = 10) +
                           s(year, bs = 'cr', k = 10) +
                           s(doy, bs = 'cc', k = 10),
                         ~ s(x, y, bs = 'ds', k = 20) +
                           s(elev_m, bs = 'ad', k = 10) +
                           s(year, bs = 'cr', k = 10) +
                           s(doy, bs = 'cc', k = 10)),
                    family = betals(),
                    knots = list(doy = c(0, 1)),
                    data = d,
                    method = 'REML',
                    control = gam.control(nthreads = 10, trace = TRUE)))
  saveRDS(m_betals, 'models/sardinia-test/betals-gam.rds')
  
  draw(m_betals, rug = FALSE)
  ggsave('figures/sardinia-test/sardinia-ndvi-betals-terms.png',
         width = 9, height = 12, units = 'in', dpi = 300, bg = 'white')
  
  summary(m_betals)
}

# using markov random fields informed with delaunay triangulation ----
# see:
# - https://groups.google.com/g/r-inla-discussion-group/c/dcrajXPn-eo
# - https://webhomes.maths.ed.ac.uk/~flindgre/posts/2018-07-22-spatially-varying-mesh-quality/
ggplot() +
  geom_sf(data = sardinia_geom) +
  geom_sf(data = locs)

# build triangular mesh
plot_mesh <- function(starting_points = NA) {
  if(is.na(starting_points)) {
    .locs <- locs
    .title <- 'Mesh with all starting points'
  } else {
    .locs <- slice_sample(locs, n = starting_points)
    .title <- paste('Mesh with', starting_points, 'random starting points')
  }
  .mesh <- fm_rcdt_2d_inla(loc = .locs,
                           boundary = fm_as_segm(sardinia_geom),
                           refine = list(max.edge = 0, min.angle = 21),
                           crs = 'EPSG:4326')
  ggplot() +
    gg(.mesh, ext.linewidth = 0.5, edge.color = 'grey',
       edge.linewidth = 0.1) +
    geom_sf(data = .locs, size = 0.75) +
    labs(x = NULL, y = NULL, title = .title)
}
plot_grid(plot_mesh(10), plot_mesh(10), plot_mesh(100), plot_mesh(),
          nrow = 1)
if(! file.exists('figures/sardinia-test/delaunay-triangulation-example.png')) {
  ggsave('figures/sardinia-test/delaunay-triangulation-example.png',
         width = 20, height = 8.68, units = 'in', dpi = 600, bg = 'white')
}

sardinia_mesh <-
  fm_rcdt_2d_inla(loc = locs,
                  boundary = fm_as_segm(sardinia_geom),
                  refine = list(max.edge = Inf, min.angle = 21),
                  crs = 'EPSG:4326')
plot(sardinia_mesh)

diagn <-
  # returns an sp object
  INLA::inla.mesh.assessment(sardinia_mesh,
                             spatial.range = 5,
                             alpha = 2,
                             dims = c(n_distinct(d$x) * 3,
                                      n_distinct(d$y) * 3)) %>%
  as.data.frame() %>%
  as_tibble() %>%
  filter(! is.na(sd)) # drop points not in a polygon

# diagnostic plots ("sd" column is tied to identifying the polygon)
diagn %>%
  filter(sd > 130, sd < 134) %>%
  ggplot() +
  coord_sf(crs = 'EPSG:4326') +
  geom_raster(aes(x, y, fill = sd)) +
  scale_fill_smoothrainbow() +
  labs(x = NULL, y = NULL) +
  theme(legend.position = 'top', legend.key.width = rel(3))

plot_grid(
  # triangle edge length is relatively similar throughout
  ggplot(diagn) +
    coord_sf(crs = 'EPSG:4326') +
    geom_raster(aes(x, y, fill = edge.len)) +
    scale_fill_viridis_c('Edge length', option = 'F', na.value = 'white') +
    labs(x = NULL, y = NULL) +
    theme(legend.position = 'top'),
  # estimate of var(pointwise) / var(continuous model); should be near 1
  ggplot(diagn) +
    coord_sf(crs = 'EPSG:4326') +
    geom_raster(aes(x, y, fill = sd.dev)) +
    scale_fill_viridis_c('Estimated variance ratio', option = 'D', na.value = 'white') +
    labs(x = NULL, y = NULL) +
    theme(legend.position = 'top'), nrow = 1)

if(! file.exists('figures/sardinia-test/delaunay-triangulation-example-diagnostics.png')) {
  ggsave('figures/sardinia-test/delaunay-triangulation-example-diagnostics.png',
         width = 8.5, height = 7.7, units = 'in', dpi = 600, bg = 'white')
}

# make list of neighboring points
sardinia_tri <- fmesher::fm_as_sfc(sardinia_mesh) %>%
  st_as_sf() %>%
  mutate(poly_id = factor(paste('poly', 1:n()))) %>%
  st_zm(drop = TRUE) # drop z axis (not used in neighbors)

ggplot() +
  geom_sf(aes(fill = poly_id), sardinia_tri) +
  scale_fill_manual(values = color('acton')(nrow(sardinia_tri))) +
  theme(legend.position = 'none')

sardinia_nbs <-
  spdep::poly2nb(pl = sardinia_tri,
                 row.names = sardinia_tri$poly_id, # so neighbors match ID
                 queen = TRUE) # 1 point in common is sufficient
names(sardinia_nbs) <- attr(sardinia_nbs, 'region.id') # names need to match names

# comparing gaussian and beta models ----
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
if(! file.exists('sardinia-ndvi-model-agreement-gaussian-beta.png')) {
  ggsave('sardinia-ndvi-model-agreement-gaussian-beta.png', plot = p_fits,
         path = 'figures/sardinia-test', width = 9, height = 6,
         units = 'in', dpi = 300, bg = 'white')
}

# plots of predicted means and squared residuals for day 1 ----
preds_1 <- d %>%
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
    labs(title = paste('Predictions for', d$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_long) +
    scale_fill_viridis_c(expression((nu-E(nu))^2), limits = c(0, NA)) +
    labs(title = paste('Squared residuals', d$date[1]), x = NULL, y = NULL),
  labels = 'AUTO', ncol = 1)

ggsave('figures/sardinia-test/sardinia-ndvi-model-agreement-gaussian-beta-predictions.png',
       height = 12, width = 8, bg = 'white')

# fit a gaussian model with markov random fields ----
ggplot(sardinia) + geom_sf(aes(fill = poly_id)) + scale_fill_bright()

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

# add cell ID and make list of neighbor cells
r_0 <- locs %>%
  st_coordinates() %>%
  data.frame() %>%
  mutate(z = 0) %>%
  rast()
plot(r_0)

sardinia_nbs <-
  adjacent(r_0, cells = cells(r_0), directions = 8, include = TRUE) %>%
  as.data.frame() %>%
  transmute(ref_cell = V1, # first column is the starting cell
            # add the 8 surrounding neighbors
            adjacent = map(1:n(), \(i) {
              .z <- c(V2[i], V3[i], V4[i], V5[i], V6[i], V7[i], V8[i], V9[i])
              
              .values <- map_lgl(.z, \(.cell_id) {
                if(is.nan(.cell_id)) {
                  return(NA)
                } else {
                  return(r_0[.cell_id]$z[1])
                }})
              
              .z <- .z[which(! is.na(.values))]
              
              if(length(.z) == 0) {
                return(0)
              } else {
                return(as.character(.z))
              }
            }))
names(sardinia_nbs$adjacent) <- sardinia_nbs$ref_cell
sardinia_nbs <- sardinia_nbs$adjacent

d <-
  mutate(d,
         cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(sardinia_nbs)))

# 
all.equal(sort(as.character(cells(r_0))), sort(names(sardinia_nbs)))

all(names(sardinia_nbs) %in% unique(d$cell_id))
all(unique(d$cell_id) %in% names(sardinia_nbs))

# all values in neighbor list are in the factor levels
all.equal(sort(levels(d$cell_id)),
          sort(as.character(unique(unlist(sardinia_nbs)))))

# all factor levels match the neighbor levels
all.equal(levels(d$cell_id), names(sardinia_nbs))

if(file.exists('models/sardinia-test/gaussian-mrf-gam.rds')) {
  m_mrf <- readRDS('models/sardinia-test/gaussian-mrf-gam.rds')
} else {
  d_1981_06_25 <- filter(d, date == '1981-06-25')
  ggplot(d_1981_06_25) +
    geom_raster(aes(x, y, fill = ndvi))
  
  m_mrf_1981_06_25 <-
    bam(ndvi_scaled ~ s(cell_id, bs = 'mrf', k = 200,
                        xt = list(nb = sardinia_nbs[unique(d_1981_06_25$cell_id)])),
        family = gaussian(),
        data = d_1981_06_25,
        method = 'fREML',
        discrete = TRUE,
        drop.unused.levels = TRUE,
        control = gam.control(trace = TRUE))
  
  plot_grid(
    plot_mrf(.model = m_mrf_1981_06_25, .term = c('(Intercept)', 's(cell_id)'),
             .newdata = d0),
    plot_mrf(.model = m_mrf_1981_06_25, .term = c('(Intercept)', 's(cell_id)'),
             .newdata = d0, .limits = c(NA, NA)),
    ggplot() +
      geom_point(aes(d0$ndvi, ndvi_to_11(fitted(m_mrf_1981_06_25))), alpha = 0.1) +
      geom_abline(intercept = 0, slope = 1, color = 'red'),
    nrow = 1)
  
  # mrf only
  m_mrf_0 <- bam(ndvi_scaled ~ s(cell_id, bs = 'mrf', xt = list(nb = nbs), k = 200),
                 family = gaussian(),
                 data = d,
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
                 data = d,
                 method = 'fREML',
                 discrete = TRUE,
                 control = gam.control(trace = TRUE))
  
  # adding a complex mrf basis informed by delaunay triangulation
  m_mrf_2 <- bam(ndvi_scaled ~
                   s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
                   s(x, y, bs = 'ds', k = 200) +
                   s(elev_m, bs = 'ad', k = 10) +
                   s(year, bs = 'cr', k = 10) +
                   s(doy, bs = 'cc', k = 10),
                 family = gaussian(),
                 knots = list(doy = c(0, 1)),
                 data = d,
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
                   data = d,
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
  tibble(ds = predict(m_gaus, newdata = d, type = 'response') %>%
           ndvi_to_11(),
         mrf = predict(m_mrf, newdata = d, type = 'response') %>%
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
preds_1_mrf <- d %>%
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
    labs(title = paste('Predictions for', d$date[1]), x = NULL, y = NULL),
  # plot squared residuals
  ggplot() +
    coord_equal() +
    facet_wrap(~ model) +
    geom_raster(aes(x, y, fill = e2), preds_1_mrf_long) +
    scale_fill_viridis_c(expression((nu-E(nu))^2), limits = c(0, NA)) +
    labs(title = paste('Squared residuals', d$date[1]), x = NULL, y = NULL),
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
                      data = d,
                      method = 'fREML',
                      discrete = TRUE,
                      control = gam.control(nthreads = 1, trace = TRUE))
  saveRDS(m_gaus_ti_ds, 'models/sardinia-test/gaus-gam-ti-ds-gam.rds')
}

draw(m_gaus_ti_ds, rug = FALSE)
ggsave('figures/sardinia-test/sardinia-ndvi-gaussian-ti-terms.png',
       width = 13.5, height = 6, units = 'in', dpi = 300, bg = 'white')

summary(m_gaus_ti_ds)
