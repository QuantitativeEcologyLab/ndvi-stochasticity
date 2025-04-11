library('dplyr')   # for data wrangling
library('sf')      # for shapefiles
library('mgcv')    # for Genralized Additive Models
library('sf')      # for simple fleature objects
library('terra')   # for rasters
library('gratia')  # for plotting GAMs
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/add_nb.R')

ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')

# set up the data for the model
if(file.exists('data/hbam-ndvi-data.rds')) {
  d <- readRDS('data/hbam-ndvi-data.rds')
} else {
  # resolution is (0.05 degrees)^2 = (3 degree minutes)^2
  # create raster of ecoregion
  r_eco <- rast('data/ecoregions/wwf-ecoregions.tif')
  
  # create raster of polygon ID
  r_poly_id <- rast('data/ecoregions/ecoregion-polygon-id.tif')
  
  # get raster of elevations
  r_elev <- rast('data/elev-raster.tif')
  
  # create raster of distance from coast
  r_dist <- rast('data/distance-from-coast-m.tif')
  
  # takes a few minutes to import
  d0 <- readRDS('data/ndvi-global-15-day-average.rds')
  
  tictoc::tic()
  #' `exactextractr::exact_extract()` fails due to bad  class for coords
  #' takes 2.3 hours on EME linux
  d <- d0 %>%
    mutate(
      year = lubridate::year(central_date),
      doy = lubridate::yday(central_date),
      poly_id = extract(r_poly_id, data.frame(x, y))[, 2], # 1 is pixel ID
      wwf_ecoregion = extract(r_eco, data.frame(x, y))[, 2],
      # distance_coast_m = extract(r_dist, data.frame(x, y)))[, 2],
      elevation_m = extract(r_elev, data.frame(x, y))[, 2]) %>%
    # convert strings to factors while separating the two hemispheres
    mutate(wwf_ecoregion = paste(wwf_ecoregion, if_else(y < 0, 'S', 'N')) %>%
             factor())
  
  d <- na.omit(d)
  nrow(d) / nrow(d0) # ~1% data loss due to NAs
  
  # some excessively low elevations for coastal pixels sea level
  # occurs mostly near islands; leaving values as they are because it helps
  # inform about the lower NDVI due to the coast
  if(FALSE) {
    plot(st_geometry(ecoregions))
    d %>%
      filter(elevation_m < -20) %>%
      points(y ~ x, ., col = 'red', pch = 19)
  }
  
  d <- mutate(d, elevation_m = if_else(elevation_m < 0, 0, elevation_m))
  
  saveRDS(d, 'data/hbam-ndvi-data.rds')
  tictoc::toc()
  
  # clean up before fitting models
  rm(d0)
  gc()
}

# create list of neighbors ----
# drop polygons with no data or neighbors
ecoregions$in_data <- ecoregions$poly_id %in% unique(d$poly_id)
filter(ecoregions, ! in_data & n_neigh == 0)

if(FALSE) {
  ggplot(ecoregions, aes(log10(area_km2), in_data)) + geom_point()
  
  ggplot(ecoregions, aes(log10(area_km2), n_neigh)) +
    facet_wrap(~ in_data) +
    geom_point()
}

ecoregions <- filter(ecoregions, in_data | n_neigh > 0)

if(file.exists('data/ecoregions/poly-nbs-global.rds')) {
  nbs <- readRDS('data/ecoregions/poly-nbs-global.rds')
} else {
  nbs <-
    spdep::poly2nb(pl = ecoregions,
                   row.names = ecoregions$poly_id, # so that names match IDs
                   queen = TRUE) # 1 point in common is sufficient
  length(nbs) == nrow(ecoregions)
  names(nbs) <- attr(nbs, 'region.id') # names do not need to match indices
  
  # add missing neighbors because of polygons that cross the edge of the map
  c('poly 8225', 'poly 8038', 'poly 7987', 'poly 8226', 'poly 7988', 'poly 8040') %in%
    ecoregions$poly_id
  
  layout(matrix(1:6, ncol = 2))
  # 7987 and 8038 are already neighbors
  # 7988 and 8040 are already neighbors
  # 7987 and 8226 do not touch
  add_nb(p1 = 'poly 8225', p2 = 'poly 8226', add = TRUE)
  add_nb(p1 = 'poly 8038', p2 = 'poly 8040', add = TRUE)
  add_nb(p1 = 'poly 7987', p2 = 'poly 7988', add = TRUE)
  add_nb(p1 = 'poly 7987', p2 = 'poly 8040', add = TRUE)
  add_nb(p1 = 'poly 7988', p2 = 'poly 8038', add = TRUE)
  add_nb('poly 8980', 'poly 8987', add = TRUE)
  layout(1)
  saveRDS(nbs, 'data/ecoregions/poly-nbs-global.rds')
}

# ensure names match the factor levels
all(attr(nbs, 'region.id') == ecoregions$poly_id)

d <- mutate(d, poly_id = factor(poly_id, levels = names(nbs)))

DATE <- format(Sys.Date(), '%Y-%m-%d')

#' *TEMPORARY VARIABLE FOR PROPORTION OF DATA*
SAMPLE <- sample(1:(nrow(d) / 3),      # only first part of the years
                 floor(nrow(d) / 1e4)) # only a subset of the rows

if(file.exists('models/global-models/hbam-mean-ndvi-DATE.rds')) {
  m <- readRDS(paste0('models/global-models/hbam-mean-ndvi-DATE.rds'))
} else {
  m_fe <- bam(
    ndvi_15_day_mean ~
      wwf_ecoregion + # to avoid intercept shrinkage
      s(poly_id, bs = 'mrf', xt = list(nb = nbs)) +
      s(doy, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cc'), k = 10) +
      s(year, wwf_ecoregion, bs = 'fs', xt = list(bs = 'cr'), k = 10) +
      ti(doy, year, wwf_ecoregion, bs = c('cc', 'cr', 're'), k = c(5, 5)) +
      s(elevation_m, bs = 'cr', k = 5),
    family = gaussian(),
    #' *TEMPORARY TEST (\/ DELETE)*
    data = d[1:(nrow(d) / 3), ],
    method = 'fREML',
    knots = list(doy = c(0.5, 366.5)),
    drop.unused.levels = TRUE,
    discrete = TRUE,
    samfrac = 0.001, # find intial guesses with a subset of the data
    nthreads = 50,
    control = gam.control(trace = TRUE))
  
  saveRDS(m_fe, paste0('models/global-models/mean-ndvi-hbam-', DATE, '.rds'))
}

#' plot each smooth separately
p_hbam_doy <- draw(m_fe, rug = FALSE,
                   select = which(grepl('s(doy', smooths(m_fe))))
ggsave(paste0('figures/hbam-mean-ndvi-fe-doy-', DATE, '.png'),
       plot = p_hbam_doy, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_y <- draw(m_fe, rug = FALSE,
                 select = which(grepl('s(year', smooths(m_fe))))
ggsave(paste0('figures/hbam-mean-ndvi-fe-year-', DATE, '.png'),
       plot = p_hbam_y, width = 24, height = 36, units = 'in', dpi = 300,
       bg = 'white')

p_hbam_elev <- draw(m_fe, rug = FALSE,
                    select = which(smooths(m_fe) == 's(elev_m)'))
ggsave(paste0('figures/hbam-mean-ndvi-fe-elev_m-', DATE, '.png'),
       plot = p_hbam_elev, width = 4, height = 6, units = 'in', dpi = 300,
       bg = 'white')

p_hbam <- draw(m, rug = FALSE, - which(smooths(m) == 's(poly_id)'))

ggsave(paste0('figures/hbam-mean-ndvi-', DATE, '.png'), plot = p_hbam,
       width = 8, height = 12, units = 'in', dpi = 600, bg = 'white')

# find pixel-level variance ----
# add fitted values and residuals to the dataset
d$mu_hat <- fitted(m)
d$e <- residuals(m)

# check pixel-level mean residuals
d %>%
  slice_sample(prop = 0.1) %>%
  group_by(long, lat) %>%
  summarize(mean_e = mean(e, na.rm = TRUE), .groups = 'drop') %>%
  pull(mean_e) %>%
  hist()

#' **IF REMOVING PIXEL-LEVEL MEAN RESIDUAL**
if(FALSE) {
  d$e <- d$e - d$mean_e
  d$mu_hat <- d$mu_hat + d$mean_e
}

# calculate pixel-level squared residual (mean e^2 is the variance)
d$e_2 <- d$e^2

saveRDS(d, 'data/hbam-var-ndvi-data.rds')
