library('dplyr')     # for data wrangling
library('tidyr')     # for data wrangling
library('sf')        # for shapefiles
library('terra')     # for rasters
library('elevatr')   # for digital elevation models
library('lubridate') # for working with dates
library('purrr')     # for functional programming
library('furrr')     # for parallelized functional programming
library('mgcv')      # for GAMs
library('ggplot2')   # for fancy plots
library('cowplot')   # for fancy plots in grids
library('gratia')    # for fancy plots of GAMs
source('functions/ndvi-palette.R')
source('functions/plot_mrf.R') # for plotting markov random field smooths
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R') # get_legend() from cowplot v. 1.1.3 fails
source('functions/nbs_from_rast.R') # gives a list of neighboring cells

# pick a northern polygon
eco <- read_sf('data/ecoregions/ecoregions-polygons.shp') %>%
  filter(WWF_REALM == 'NA') %>%
  st_transform(crs(rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc')))
plot(st_geometry(eco), col = 'grey')
plot(st_geometry(eco)[18, ], col = 'red', add = TRUE)

arctic <- eco %>%
  slice(18) %>%
  st_cast('POLYGON', warn = FALSE) %>%
  st_geometry() %>%
  st_as_sf() %>%
  slice(100)
plot(st_geometry(eco), col = 'grey')
plot(arctic, col = 'red', add = TRUE)

# there are some oddly large NDVI values in winter
list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
           pattern = '.nc', full.names = TRUE)[seq(1, 365, by = 10)] %>%
  future_map(\(fn) {
    r <- rast(fn, lyr = 'NDVI')
    r %>%
      crop(., st_transform(arctic, crs(.))) %>%
      mask(st_transform(arctic, crs(.))) %>%
      as.data.frame(xy = TRUE, na.rm = TRUE) %>%
      mutate(date = as.Date(unique(time(r))))
  }, .progress = TRUE) %>%
  bind_rows() %>%
  mutate(doy = yday(date)) %>%
  ggplot(aes(doy, NDVI)) +
  geom_point(alpha = 0.1) +
  geom_smooth(formula = y ~ s(x, k = 5, bs = 'cc'), method = 'gam',
              method.args = list(knots = list(doy = c(0.5, 366.5))),
              fullrange = TRUE) +
  xlim(c(0, 366))

# import ndvi data ----
if(file.exists('data/arctic-test/arctic-ndvi.rds')) {
  d <- readRDS('data/arctic-test/arctic-ndvi.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  d <-
    list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI')
      r %>%
        crop(., st_transform(arctic, crs(.))) %>%
        mask(st_transform(arctic, crs(.))) %>%
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI)
  plan(sequential)
  
  # the period with NDVI > 0.1 is fairly narrow
  ggplot(filter(d, date %in% unique(date)[1:21]), aes(date, ndvi)) +
    geom_point() +
    geom_smooth(formula = y ~ s(x, k = 5))
  
  # plot the first few NDVI rasters
  p_d <- filter(d, date %in% unique(date)[1:21]) %>%
    ggplot() +
    facet_wrap(~ date, nrow = 3) +
    geom_raster(aes(x, y, fill = ndvi)) +
    geom_sf(data = arctic, fill = 'transparent') +
    scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
    ylab(NULL) +
    scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
  
  # add altitude (get_elev_points results in a json error)
  elevs <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
    crop(st_buffer(arctic, 1e4)) %>% # crop to area
    rast() # convert rasterLayer to SpatRaster
  
  plot(elevs)
  plot(arctic, add = TRUE, col = 'transparent', lwd = 2)
  
  d <- mutate(d,
              elev_m = extract(elevs, select(d, x, y))[, 2],
              year = year(date),
              doy = yday(date))
  range(d$elev_m)
  quantile(d$elev_m, c(0.1, 0.01, 0.001))
  
  d <- mutate(d, elev_m = if_else(elev_m < 0, 0, elev_m))
  
  d <- filter(d, ! is.na(cell_id)) # drop rows with no cell id
  
  saveRDS(d, 'data/arctic-test/arctic-ndvi.rds')
}

summary(d) # max(ndvi) = 1!

ggplot(slice_sample(d, n = 1e5), aes(doy, ndvi)) +
  geom_point(alpha = 0.1) +
  geom_smooth(formula = y ~ s(x, k = 10, bs = 'cc'), method = 'gam',
              method.args = list(knots = list(doy = c(0.5, 366.5))),
              fullrange = TRUE) +
  xlim(c(0, 366))

bam(ndvi ~ s(doy, bs = 'cc'), data = slice_sample(d, n = 1e5),
    method = 'REML', knots = list(doy = c(0.5, 366.5))) %>%
  plot(scheme = 1, xlim = c(0, 366), n = 400, resid = TRUE)

# add cell ID and make list of neighbor cells for all coordinates ----
# unique coordinates inside arctic
locs <- d %>%
  select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  filter(., st_as_sf(., coords = c('x', 'y')) %>%
           st_set_crs('EPSG:4326') %>%
           st_intersects(st_transform(arctic, 'EPSG:4326'),
                         sparse = TRUE) %>%
           map_lgl(\(x) length(x) > 0))
nrow(locs) # all rasters use same coords

p_locs <-
  ggplot() +
  geom_sf(data = arctic) +
  geom_sf(data = locs)
p_locs

# make a raster of all locations
r_0 <- locs %>%
  st_coordinates() %>%
  data.frame() %>%
  mutate(z = 0) %>%
  rast()

p_locs +
  geom_raster(aes(x, y), as.data.frame(r_0, xy = TRUE), fill = '#FF000030') +
  labs(x = NULL, y = NULL)

nbs <- nbs_from_rast(r_0)

d <-
  mutate(d,
         cell_id = cellFromXY(r_0, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs)))

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0))), sort(names(nbs)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
all.equal(sort(levels(d$cell_id)),
          sort(as.character(unique(unlist(nbs)))))

# fit a spatially explicit test model with a gaussian family ----
if(all(file.exists(c('models/arctic-test/gaussian-gam-ds.rds',
                     'models/arctic-test/gaussian-gam-mrf.rds')))) {
  m_gaus_ds <- readRDS('models/arctic-test/gaussian-gam-ds.rds')
  m_gaus_mrf <- readRDS('models/arctic-test/gaussian-gam-mrf.rds')
} else {
  m_gaus_ds <- bam(
    ndvi ~
      s(x, y, bs = 'ds', k = 200) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0.5, 366.5)),
    data = d,
    method = 'fREML',
    discrete = TRUE,
    nthreads = 10)
  
  m_gaus_mrf <- bam(
    ndvi ~
      s(cell_id, bs = 'mrf', k = 200, xt = list(nb = nbs)) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0.5, 366.5)),
    data = d,
    method = 'fREML',
    discrete = TRUE)
  
  # deviance explained and complexity of the spatial terms are similar
  summary(m_gaus_ds)
  summary(m_gaus_mrf)
  
  draw(m_gaus_ds, rug = FALSE, dist = 0.07)
  ggsave('figures/arctic-test/arctic-ndvi-gaussian-ds-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  plot_mrf(.model = m_gaus_mrf, .newdata = d, .full_model = TRUE)
  ggsave('figures/arctic-test/arctic-ndvi-gaussian-mrf-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  
  saveRDS(m_gaus_ds, 'models/arctic-test/gaussian-gam-ds.rds')
  saveRDS(m_gaus_mrf, 'models/arctic-test/gaussian-gam-mrf.rds')
}

# testing data aggregation ----
s_res <- 4 # spatial resolution
t_res <- 4 # temporal resolution

# import and aggregate the data ----
if(file.exists('data/arctic-test/arctic-ndvi-t-4-s-4-aggr.rds')) {
  d_aggr <- readRDS('data/arctic-test/arctic-ndvi-t-4-s-4-aggr.rds')
} else {
  if(.Platform$OS.type != 'unix')
    stop('AVHRR/VIIRS rasters are on the H: Drive, and you may want to use multiple cores.')
  future::availableCores(logical = FALSE)
  plan(multisession, workers = min(60, availableCores(logical = FALSE) - 2))
  d_aggr <-
    list.files(path = 'data/avhrr-viirs-ndvi/raster-files/', #'/home/shared/NOAA_Files/',
               pattern = '.nc',
               full.names = TRUE) %>%
    future_map(\(fn) {
      r <- rast(fn, lyr = 'NDVI') # to extract time below
      r %>%
        crop(., st_transform(arctic, crs(.))) %>%
        terra::aggregate(s_res, na.rm = TRUE) %>%
        mask(st_transform(arctic, crs(.))) %>% # mask after aggregating
        as.data.frame(xy = TRUE, na.rm = TRUE) %>%
        mutate(date = as.Date(unique(time(r))))
    }, .progress = TRUE) %>%
    bind_rows() %>%
    as_tibble() %>%
    rename(ndvi = NDVI) %>%
    # aggregate temporally
    mutate(julian = julian(date),
           central_date = as.Date(julian - (julian %% t_res) + 2)) %>%
    group_by(central_date, x, y) %>%
    summarize(doy = yday(central_date),
              year = year(central_date),
              ndvi = mean(ndvi, na.rm = TRUE)) %>%
    ungroup()
  plan(sequential)
  
  # plot the first few NDVI rasters
  p_d_aggr <-
    filter(d_aggr, central_date %in% unique(as.character(central_date))[1:21]) %>%
    ggplot() +
    facet_wrap(~ as.character(central_date), nrow = 3) +
    geom_raster(aes(x, y, fill = ndvi)) +
    geom_sf(data = arctic, fill = 'transparent') +
    scale_x_continuous(NULL, breaks = c(8.5, 9.5)) +
    ylab(NULL) +
    scale_fill_gradientn('NDVI', colors = ndvi_pal, limits = c(-1, 1))
  p_d_aggr
  
  # add altitude (get_elev_points results in a json error)
  elevs_aggr <- d %>%
    filter(date == first(date)) %>%
    select(x, y) %>%
    mutate(z = 1) %>%
    rast(crs = 'EPSG:4326') %>%
    aggregate(4, na.rm = TRUE) %>%
    get_elev_raster(z = 2) %>% # nearest finer res than 0.20x0.20
    crop(st_buffer(arctic, 1e4))
  
  plot(elevs_aggr)
  plot(arctic, add = TRUE, col = 'transparent')
  
  d_aggr <- mutate(d_aggr,
                   elev_m = extract(elevs_aggr, select(d_aggr, x, y)),
                   year = year(central_date),
                   doy = yday(central_date))
  
  range(d_aggr$elev_m)
  quantile(d_aggr$elev_m, c(0.1, 0.01, 0.001))
  
  saveRDS(d_aggr, paste0('data/arctic-test/arctic-ndvi-t-',
                         t_res, '-s-', s_res, '-aggr.rds'))
}

# ndvi still has high max values
summary(d_aggr)

# plot a comparison of the first 21 rasters
plot_grid(
  get_legend(p_d +
               theme(legend.position = 'top', legend.key.width = rel(2))),
  plot_grid(p_d + theme(legend.position = 'none'),
            p_d_aggr + theme(legend.position = 'none'),
            labels = 'AUTO'),
  rel_heights = c(1, 15), ncol = 1)

ggsave('figures/arctic-test/example-rasters-aggregation.png',
       width = 20, height = 8.5, units = 'in', dpi = 600, bg = 'white')

# fit a model to show the changes in predictions ----
# add cell ID and make list of neighbor cells for all coordinates ----
# unique coordinates inside arctic
locs_aggr <- d_aggr %>%
  select(x, y) %>%
  group_by(x, y) %>%
  slice(1) %>%
  st_as_sf(coords = c('x', 'y')) %>%
  st_set_crs('EPSG:4326') %>%
  filter(., st_as_sf(., coords = c('x', 'y')) %>%
           st_set_crs('EPSG:4326') %>%
           st_transform(crs(arctic)) %>%
           st_intersects(arctic, sparse = TRUE) %>%
           map_lgl(\(x) length(x) > 0))
nrow(locs_aggr) # all rasters use same coords

p_locs_aggr <-
  ggplot() +
  geom_sf(data = arctic) +
  geom_sf(data = locs_aggr)
p_locs_aggr

# make a raster of all locations
r_0_aggr <- locs_aggr %>%
  st_coordinates() %>%
  data.frame() %>%
  mutate(z = 0) %>%
  rast()

p_locs_aggr +
  geom_raster(aes(x, y), as.data.frame(r_0_aggr, xy = TRUE),
              fill = '#FF000030') +
  labs(x = NULL, y = NULL)

nbs_aggr <- nbs_from_rast(r_0_aggr)

d_aggr <-
  mutate(d_aggr,
         cell_id = cellFromXY(r_0_aggr, xy = as.matrix(tibble(x, y))) %>%
           factor(levels = names(nbs_aggr))) %>%
  filter(! is.na(cell_id))
mean(is.na(d_aggr$cell_id)) # some cell centers fall outside of the polygon

n_distinct(d_aggr$cell_id)
length(names(nbs_aggr))

# all cell names match the neighbor list names
# there cannot be any list elements with names that are not in the dataset
all.equal(sort(as.character(cells(r_0_aggr))), sort(names(nbs_aggr)))

# all values in neighbor list are in the factor levels
# this is not crucial, as the model will predict for neighbors with data
all.equal(sort(levels(d_aggr$cell_id)),
          sort(as.character(unique(unlist(nbs_aggr)))))

if(file.exists('models/arctic-test/gaussian-gam-ds.rds')) {
  m_gaus_mrf_aggr <- readRDS('models/arctic-test/gaussian-gam-mrf-aggr.rds')
} else {
  m_gaus_mrf_aggr <- bam(
    ndvi ~
      s(cell_id, bs = 'mrf', k = 14, xt = list(nb = nbs_aggr)) +
      s(elev_m, bs = 'cr', k = 5) +
      s(year, bs = 'cr', k = 10) +
      s(doy, bs = 'cc', k = 10),
    family = gaussian(),
    knots = list(doy = c(0, 1)),
    data = d_aggr,
    method = 'fREML',
    discrete = TRUE)
  plot_mrf(.model = m_gaus_mrf_aggr, .newdata = d_aggr, .full_model = TRUE)
  ggsave('figures/arctic-test/arctic-ndvi-gaussian-mrf-aggr-terms.png',
         width = 9, height = 6, units = 'in', dpi = 300, bg = 'white')
  saveRDS(m_gaus_mrf_aggr, 'models/arctic-test/gaussian-gam-mrf-aggr.rds')
}

# make a figure comparing ds gam, mrf gam, and mrf_aggr gam ----
elevs <- d %>%
  filter(date == first(date)) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  get_elev_raster(z = 4) %>% # nearest finer res than 0.05x0.05
  crop(st_buffer(arctic, 1e4)) # crop to area near arctic

elevs_aggr <- d %>%
  filter(date == first(date)) %>%
  select(x, y) %>%
  mutate(z = 1) %>%
  rast(crs = 'EPSG:4326') %>%
  aggregate(4, na.rm = TRUE) %>%
  get_elev_raster(z = 2) %>% # nearest finer res than 0.20x0.20
  crop(st_buffer(arctic, 1e4))

gratia::smooths(m_gaus_ds)
gratia::smooths(m_gaus_mrf)

# get model predictions ----
get_preds <- function(nd, space = TRUE) {
  if(space) {
    preds <- nd %>%
      mutate(.,
             ds_mu =
               predict(object = m_gaus_ds, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(x,y)', 's(elev_m)')),
             mrf_mu =
               rename(., cell_id = cell_id_fine) %>%
               predict(object = m_gaus_mrf, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')),
             mrfa_mu =
               rename(., cell_id = cell_id_aggr) %>%
               predict(object = m_gaus_mrf_aggr, newdata = .,
                       type = 'response', se.fit = FALSE,
                       terms = c('(Intercept)', 's(cell_id)', 's(elev_m)')))
  } else {
    preds <- nd %>% # new data
      mutate(
        ds_mu = predict(object = m_gaus_ds, newdata = .,
                        type = 'response', se.fit = FALSE,
                        terms = c('(Intercept)', 's(doy)')),
        mrf_mu = predict(object = m_gaus_mrf, newdata = .,
                         type = 'response', se.fit = FALSE,
                         terms = c('(Intercept)', 's(doy)')),
        mrfa_mu = predict(object = m_gaus_mrf_aggr, newdata =.,
                          type = 'response', se.fit = FALSE,
                          terms = c('(Intercept)', 's(doy)')))
  }
  
  preds <- preds %>%
    mutate(diff_mu = mrfa_mu - mrf_mu) %>%
    select(x, y, elev_m, doy, ds_mu:diff_mu)
  
  if(space) {
    # get temporally static maps of estimated variance
    get_s2 <- function(.model) {
      .m <- if(.model == 'Duchon spline') {
        .m <- m_gaus_ds
      } else if(.model == 'MRF') {
        .m <- m_gaus_mrf
      } else if(.model == 'MRF (aggregated)') {
        .m <- m_gaus_mrf_aggr
      }
      
      # need to select the dataset based on the model but keep (x,y) coords
      .d <- if(.model == 'Duchon spline' | .model == 'MRF') {
        .d <- d
      } else if(.model == 'MRF (aggregated)') {
        .d <- na.omit(d_aggr)
      }
      
      .d %>%
        transmute(x, y, e2 = resid(.m)^2) %>%
        group_by(x, y) %>%
        summarise(s2 = mean(e2), .groups = 'drop') %>%
        rast() %>%
        `crs<-`('EPSG:4326') %>%
        project(elevs, res = res(elevs)) %>%
        extract(., select(as.data.frame(elevs, xy = TRUE), 1:2)) %>%
        select(! ID) %>%
        bind_cols(select(as.data.frame(elevs, xy = TRUE), 1:2), .) %>%
        mutate(model = .model) %>%
        filter(! is.na(s2)) %>%
        as_tibble() %>%
        mutate(model = .model) %>%
        return()
    }
    
    s2 <- bind_rows(map(c('Duchon spline', 'MRF', 'MRF (aggregated)'), get_s2)) %>%
      pivot_wider(names_from = model, values_from = s2) %>%
      mutate(diff = `MRF (aggregated)` - `MRF`) %>%
      pivot_longer(`Duchon spline`:diff, names_to = 'model', values_to = 'value')
    
    preds <-
      bind_rows(pivot_longer(preds, ds_mu:diff_mu,
                             names_to = c('model', 'param'),
                             values_to = 'value', names_sep = '_'),
                mutate(s2, param = 's2'))
  } else {
    preds <- 
      left_join(
        preds,
        tibble(
          doy = unique(d$doy),
          ds_s2 = d %>%
            mutate(e2 = resid(m_gaus_ds)^2) %>%
            group_by(doy) %>%
            summarize(s2 = mean(e2, na.rm = TRUE)) %>%
            pull(s2),
          mrf_s2 = d %>%
            mutate(e2 = resid(m_gaus_mrf)^2) %>%
            group_by(doy) %>%
            summarize(s2 = mean(e2, na.rm = TRUE)) %>%
            pull(s2)),
        by = 'doy') %>%
      left_join(
        tibble(doy = unique(d_aggr$doy),
               mrfa_s2 = d_aggr %>%
                 na.omit() %>%
                 mutate(e2 = resid(m_gaus_mrf_aggr)^2) %>%
                 group_by(doy) %>%
                 summarize(s2 = mean(e2, na.rm = TRUE)) %>%
                 pull(s2)),
        by = 'doy') %>%
      mutate(diff_s2 = mrfa_s2 - mrf_s2) %>%
      select(doy, ds_mu:diff_s2) %>%
      pivot_longer(ds_mu:diff_s2,
                   names_to = c('model', 'param'), values_to = 'value',
                   names_sep = '_')
  }
  preds %>%
    mutate(model = case_when(model == 'ds' ~ 'Duchon spline',
                             model == 'mrf' ~ 'MRF',
                             model == 'mrfa' ~ 'MRF (aggregated)',
                             model == 'diff' ~ 'diff',
                             TRUE ~ model)) %>%
    return()
}

# spatial predictions
preds_comp_s <-
  elevs %>%
  mask(arctic) %>%
  as.data.frame(xy = TRUE) %>%
  rename(elev_m = 3) %>%
  filter(! is.na(elev_m)) %>%
  mutate(elev_m = if_else(elev_m < 0, 0, elev_m)) %>%
  mutate(cell_id_fine = factor(cells(r_0, vect(tibble(x, y),
                                               geom = c('x', 'y')))[, 2],
                               levels = levels(d$cell_id)),
         cell_id_aggr = factor(cells(r_0_aggr,
                                     vect(tibble(x, y),
                                          geom = c('x', 'y')))[, 2],
                               levels = levels(d_aggr$cell_id)),
         year = 0, doy = 0) %>%
  get_preds(space = TRUE) # add model predictions

# temporal doy predictions
preds_comp_t <- tibble(doy = 1:366,
                       cell_id = intersect(d$cell_id, d_aggr$cell_id)[1],
                       x = 0, y = 0, elev_m = 0, year = 0) %>%
  get_preds(space = FALSE)

# make the final figure ----
p_comp <-
  plot_grid(
    ncol = 2, labels = 'AUTO', rel_widths = c(3, 1.5),
    rel_heights = c(3, 3, 2, 2, 2),
    # row 1: map of mean NDVI
    ggplot(filter(preds_comp_s, param == 'mu', model != 'diff')) +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = arctic, fill = 'transparent', color = 'black') +
      scale_fill_viridis_c('NDVI', option = 'A') +
      labs(x = NULL, y = NULL),
    ggplot(filter(preds_comp_s, param == 'mu', model == 'diff')) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = arctic, fill = 'transparent', color = 'black') +
      scale_fill_distiller(
        expression(atop(bold('Difference in'), bold('mean NDVI'))),
        type = 'div', palette = 5,
        limits = max(abs(filter(preds_comp_s, model == 'diff',
                                param == 'mu')$value)) * c(-1, 1)) +
      labs(x = NULL, y = NULL),
    # row 2: map of variance in NDVI
    filter(preds_comp_s, param == 's2', model != 'diff') %>%
      group_by(model) %>%
      mutate(value = value / max(value, na.rm = TRUE)) %>%
      ggplot() +
      facet_grid(. ~ model) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = arctic, fill = 'transparent', color = 'black') +
      scale_fill_viridis_c(expression(bold(s^'2'))) +
      labs(x = NULL, y = NULL),
    ggplot(filter(preds_comp_s, param == 's2', model == 'diff')) +
      geom_raster(aes(x, y, fill = value)) +
      geom_sf(data = arctic, fill = 'transparent', color = 'black') +
      scale_fill_distiller(
        expression(bold(atop('Difference in', 'estimated s'^'2'))),
        type = 'div', palette = 4) +
      labs(x = NULL, y = NULL),
    # row 3: mean NDVI over day of year
    ggplot(filter(preds_comp_t, param == 'mu', model != 'diff')) +
      facet_grid(. ~ model) +
      geom_line(aes(doy, value)) +
      labs(x = 'Day of year', y = 'Mean NDVI'),
    ggplot(filter(preds_comp_t, param == 'mu', model == 'diff')) +
      geom_line(aes(doy, value)) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      labs(x = 'Day of year', y = 'Difference in mean NDVI'),
    # row 4: variance in NDVI over day of year
    ggplot(filter(preds_comp_t, param == 's2', model != 'diff'),
           aes(doy, value)) +
      facet_grid(. ~ model) +
      geom_point(alpha = 0.3) +
      geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
                  method.args = list(knots = list(x = c(0.5, 366.5)))) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      labs(x = 'Day of year',
           y = expression(bold(paste('Daily mean e'^'2')))),
    ggplot(filter(preds_comp_t, param == 's2', model == 'diff'),
           aes(doy, value)) +
      geom_point(alpha = 0.3) +
      geom_smooth(formula = y ~ s(x, bs = 'cc'), method = 'gam',
                  method.args = list(knots = list(x = c(0.5, 366.5)))) +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      geom_rug(aes(x = x), tibble(x = unique(d$doy)), alpha = 0.1,
               inherit.aes = FALSE) +
      labs(x = 'Day of year',
           y = expression(bold(paste('Difference in daily mean e'^'2')))),
    # row 5: NDVI over elevation
    ggplot(filter(preds_comp_s, param == 'mu', model != 'diff') %>%
             filter(! is.na(value))) +
      facet_grid(. ~ model) +
      geom_point(aes(elev_m, value), alpha = 0.1) +
      geom_rug(aes(x = elev_m), alpha = 0.1) +
      labs(x = 'Elevation (m)', y = 'Mean NDVI'),
    ggplot(filter(preds_comp_s, param == 'mu', model == 'diff') %>%
             filter(! is.na(value))) +
      geom_point(aes(elev_m, value), alpha = 0.1) +
      geom_smooth(aes(elev_m, value), formula = y ~ s(x, k = 5),
                  method = 'gam') +
      geom_hline(yintercept = 0, color = 'grey', lty = 'dashed') +
      geom_rug(aes(x = elev_m), alpha = 0.01) +
      labs(x = 'Elevation (m)', y = 'Difference in mean NDVI'))

ggsave('figures/arctic-test/model-comparisons.png',
       width = 12.5, height = 17, units = 'in', dpi = 300, bg = 'white')
