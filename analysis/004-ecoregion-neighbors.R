library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('purrr')   # for data wrangling
library('sf')      # for spatial data
library('terra')   # for rasters
library('spdep')   # for finding neighbors
library('ggplot2') # for fancy plots
source('functions/add_nb.R')

# get shapefile of each continent
ecoregions <- st_read('data/world-ecosystems/data/commondata/data0/tnc_terr_ecoregions.shp') %>%
  filter(WWF_REALM2 != 'Antarctic') %>% # drop Antarctica
  #' **\/ for testing **
  # filter(WWF_REALM == 'PA') %>%
  # st_make_valid() %>%
  # st_intersection(filter(world, name_long == 'Italy') %>%
  #                   st_buffer(10e3) %>%
  #                   st_bbox() %>%
  #                   st_as_sfc() %>%
  #                   st_as_sf() %>%
  #                   st_set_crs('EPSG:4326')) %>%
  #' **/\ for testing **
  mutate(group = case_when(WWF_REALM %in% c('NA', 'NT') ~ 'americas',
                           WWF_REALM %in% c('AT', 'PA', 'IM') ~ 'africa+eurasia',
                           WWF_REALM %in% c('AA', 'OC') ~ 'islands',
                           .default = 'unassigned')) %>%
  select(group, ECO_NAME, ECO_NUM, WWF_MHTNAM) %>%
  #' casting to polygon because some ecoregions are multipolygons, and we
  #' we need a list of individual polygons for MRF smooths
  #' turning off spherical geometry to prevent `st_make_valid()` from making
  #' russia span the whole northern hemisphere by inverting it
  st_cast('POLYGON') %>% # assigns attributes to sub-geometries by default
  st_make_valid() %>% # to drop duplicate vertices
  mutate(poly_id = paste('poly', 1:n()),
         area_km2 = as.numeric(st_area(.)) / 1e6) %>%
  filter(WWF_MHTNAM != 'Inland Water') # dropping 17 inland water polygons

# drop points based on plot(n_data ~ area)
# drop points with areas < (5 km)^2 * 1, 4, or 9 cells
plot(ecdf(log10(ecoregions$area_km2)))
abline(v = log10(1), col = 'black', lty = 'dashed')
abline(v = log10(5^2 * 1), col = '#004488') # ~ 1 cell
abline(v = log10(5^2 * 4), col = '#DDAA33') # ~ 2*2 cells
abline(v = log10(5^2 * 9), col = '#BB5566') # ~ 3*3 cells
nrow(ecoregions)

# some exploratory plots and checks
if(FALSE) {
  # find total area by group to check proprotion of NDVI data per group
  ecoregions %>%
    st_drop_geometry() %>%
    group_by(group) %>%
    summarize(total = sum(area_km2)) %>%
    mutate(total = total / sum(total))
  
  # each ecoregion has a unique name, but many ecoregions are multipolygons
  unique(ecoregions$ECO_NAME)[1:10]
  plot(st_geometry(ecoregions), col = factor(ecoregions$ECO_NAME),
       axes = TRUE)
  
  # grouping is ok, but some polygons aren't continuous with neighbors
  ggplot(ecoregions) +
    facet_wrap(~ group, ncol = 1) +
    geom_sf(fill = 'red3', color = 'black') +
    theme_classic()
}

# find neighbors (warnings are ok: islands don't have neighbors)
#' `snap` has no effect with spherical geometries, which the ecoregions use
nbs <- poly2nb(pl = ecoregions,
               row.names = ecoregions$poly_id,# so that neighbors match ID
               queen = TRUE) # 1 point in common is sufficient

#' `nbs` is a list with `c(0)` if no neighbors
str(nbs[100:106]) # ith element is the integer index of the neighbors
length(nbs) == nrow(ecoregions)

# add names corresponding to each polygon
all(attr(nbs, 'region.id') == ecoregions$poly_id)
names(nbs) # no names, currently
names(nbs) <- attr(nbs, 'region.id') # names do not need to match indices
names(nbs)

# add number of neighbors for each polygon
ecoregions$n_neigh <- map_int(nbs, function(.nb) {
  if(length(.nb) == 1 & all(.nb == 0)) {
    return(0)
  } else {
    return(length(.nb))
  }
})

# drop small islands with no neighbors ----
# 0.05 deg is approx 5.5 km at the equator and 1.43 km^2 at 75 deg N
# see https://www.opendem.info/arc2meters.html
mean(ecoregions$n_neigh == 0) # some will have no nbs but enough data
mean(ecoregions$n_neigh == 0 & ecoregions$area_km2 < 5^2 * 1) # 1 cell
mean(ecoregions$n_neigh == 0 & ecoregions$area_km2 < 5^2 * 4) # 2*2 cells
mean(ecoregions$n_neigh == 0 & ecoregions$area_km2 < 5^2 * 9) # 3*3 cells

MIN_AREA <- 5^2 * 1 # (25 km)^2
nrow(ecoregions)
ecoregions <- ecoregions %>%
  filter(! (ecoregions$n_neigh == 0 & ecoregions$area_km2 < MIN_AREA))
nrow(ecoregions)

nbs <- poly2nb(pl = ecoregions,
               row.names = ecoregions$poly_id, # so that neighbors match ID
               queen = TRUE) # 1 point in common is sufficient

# update the numbers and names after dropping the islands
length(nbs)
length(nbs) == nrow(ecoregions)

all(attr(nbs, 'region.id') == ecoregions$poly_id)
names(nbs) <- attr(nbs, 'region.id') # re-add names

ecoregions$n_neigh <- map_int(nbs, function(.nb) {
  if(length(.nb) == 1 & all(.nb == 0)) {
    return(0)
  } else {
    return(length(.nb))
  }
})

# no excessively small islands
sum(ecoregions$n_neigh == 0 & ecoregions$area_km2 < MIN_AREA)

# some polys at long = 180 don't count as neighbors ----
ha_proj <- 'ESRI:102007' # hawai'i albers

polys_180 <- ecoregions %>%
  mutate(min_x = map_dbl(geometry, \(.g) st_bbox(.g)['xmin']),
         max_x = map_dbl(geometry, \(.g) st_bbox(.g)['xmax']),
         n_neigh = map_int(nbs, length)) %>%
  filter(min_x < -179 | max_x > 179) %>%
  st_as_sf()
polys_180$n_neigh
polys_180$area_km2

polys_180 %>%
  st_transform(ha_proj) %>%
  ggplot() +
  geom_sf(aes(fill = n_neigh), color = 'black') +
  scale_fill_viridis_c(direction = -1) +
  theme_minimal()

pal <- c('#004488', '#DDAA33', '#BB5566', 'grey90', 'black', 'grey50')

# add neighbors manually
ru <- polys_180 %>%
  st_as_sf() %>%
  st_transform(ha_proj) %>%
  filter(group == 'africa+eurasia')
ru_geom <- st_geometry(ru)

layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE))
plot(ru_geom, col = pal, main = 'Russian polygons that cross long = 180')
plot(ru_geom[1:4, ], col = pal)
plot(ru_geom[5:6, ], col = pal[5:6])
layout(1)

layout(matrix(1:6, ncol = 2))
plot(ru_geom, col = pal, main = 'Russian polygons that cross long = 180')
# 7987 and 8038 area already neighbors
# 7988 and 8040 are already neighbors
# 7987 and 8226 do not touch
add_nb(p1 = 'poly 8225', p2 = 'poly 8226', add = TRUE)
add_nb(p1 = 'poly 8038', p2 = 'poly 8040', add = TRUE)
add_nb(p1 = 'poly 7987', p2 = 'poly 7988', add = TRUE)
add_nb(p1 = 'poly 7987', p2 = 'poly 8040', add = TRUE)
add_nb(p1 = 'poly 7988', p2 = 'poly 8038', add = TRUE)
layout(1)

# quick check to make sure the IDs match
if(FALSE) {
  nbs[['poly 8226']]
  nbs[['poly 8225']]
  which(names(nbs) == 'poly 8226')
  which(names(nbs) == 'poly 8225')
}

# repeat for one island in the souther pacific
isl <- polys_180 %>%
  st_as_sf() %>%
  st_transform(ha_proj) %>%
  filter(group != 'africa+eurasia')
unique(isl$group)
isl_geom <- st_geometry(isl)

plot(isl_geom, col = 1:nrow(isl), main = 'Pacific islands that cross long = 180')
ggplot(isl) +
  geom_sf(aes(fill = poly_id)) +
  khroma::scale_fill_discreterainbow(name = 'Polygon ID')

add_nb('poly 8980', 'poly 8987', add = TRUE)

# update number of neighbors after adding some manually ----
sum(ecoregions$n_neigh)
ecoregions$n_neigh <- map_int(nbs, function(.nb) {
  if(length(.nb) == 1 & all(.nb == 0)) {
    return(0)
  } else {
    return(length(.nb))
  }
})
sum(ecoregions$n_neigh) # (5 + 1) * 2 additions

ggplot(ecoregions) +
  geom_sf(aes(fill = n_neigh)) +
  khroma::scale_fill_smoothrainbow(name = 'Number of neighboring polygons') +
  theme(legend.position = 'top')

# save the final files ----
st_write(ecoregions, 'data/ecoregions/ecoregions-polygons.shp')

# run a test ----
if(FALSE) {
  rstudioapi::restartSession()
  library('dplyr')   # for data wrangling
  library('tidyr')   # for data wrangling
  library('sf')      # for spatial data
  library('terra')   # for rasters
  library('ggplot2') # for fancy plots
  library('mgcv')    # for GAMs
  
  ecoregions <- st_read('data/ecoregions/ecoregions-polygons.shp')
  
  r0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc', lyr = 'NDVI') %>%
    aggregate(2) %>%
    `values<-`(1) %>%
    mask(vect(ecoregions))
  plot(r0)
  
  ecoregions$n_cells <- r0 %>%
    extract(vect(ecoregions), fun = \(x) length(x)) %>%
    pull(NDVI)
  
  d_sample <- bind_rows(readRDS('data/ndvi-global-15-day-average-2020-only.rds'),
                        readRDS('data/ndvi-global-15-day-average-2021-only.rds')) %>%
    #' use `exact_extract()` for more data?
    mutate(doy = lubridate::yday(central_date),
           poly_id = extract(r_poly_id, data.frame(x, y))[, 2], # 1 is pixel ID
           wwf_ecoregion = extract(r_eco, data.frame(x, y))[, 2],
           distance_coast_m = extract(r_dist, data.frame(x, y))[, 2],
           elevation_m = extract(r_elev, data.frame(x, y))[, 2]) %>%
    # convert strings to factors
    mutate(poly_id = factor(poly_id, levels = names(nbs)),
           wwf_ecoregion = factor(wwf_ecoregion)) %>%
    na.omit()
  
  RES <- 0.1 # AVHRR NDVI rasters have a resolution of (0.05 degrees)^2
  r <- expand_grid(x = seq(-180, 180, by = RES),
                   y = seq(-89, 89, by = RES)) %>%
    mutate(z = 1L) %>%
    rast() %>%
    `crs<-`('EPSG:4326') %>%
    crop(ecoregions) %>%
    mask(ecoregions, touches = FALSE)
  plot(r)
  
  d <- r %>%
    as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    mutate(.,
           poly_i = st_as_sf(., coords = c('x', 'y')) %>% # create index
             st_set_crs('EPSG:4326') %>%
             st_intersects(ecoregions, sparse = TRUE),
           poly_i = purrr::map_int(poly_i, function(.p) {
             # some raster cells are assigned to 0 or > 2 polygons (rarely)
             if(length(.p) == 0) {
               return(NA_integer_)
             } else if(length(.p) == 1) {
               return(as.integer(.p))
             } else {
               warning('Selected the first of ', length(.p), ' polygons.')
               return(as.integer(.p)[1])
             }
           }),
           eco_name = factor(ecoregions$ECO_NAME[poly_i]),
           poly_id = ecoregions$poly_id[poly_i],
           poly_id = factor(poly_id, levels = names(nbs)),
           z = rnorm(n = n(), mean = as.numeric(eco_name) / 500)) %>%
    select(! poly_i) %>% # remove polygon index
    as_tibble() %>%
    filter(! is.na(poly_id)) # ensure all points are assigned to a polygon
  unique(dplyr::last_dplyr_warnings())
  d
  
  # polygon names match must
  all(levels(d$poly_id) %in% names(nbs))
  all(names(nbs) %in% levels(d$poly_id))
  
  # ensure names match
  all.equal(sort(names(nbs)), sort(levels(d$poly_id)))
  
  if(FALSE) {
    # plot points with ecoregions
    ggplot() +
      geom_sf(aes(fill = ECO_NAME), data = ecoregions) +
      geom_sf(data = d %>%
                slice_sample(n = 1e4) %>%
                st_as_sf(coords = c('x', 'y')) %>%
                st_set_crs('EPSG:4326'), size = 0.05) +
      theme(legend.position = 'none')
    
    # plot points with poligon IDs
    ggplot() +
      geom_sf(aes(fill = poly_id), data = ecoregions) +
      geom_sf(data = d %>%
                slice_sample(n = 1e4) %>%
                st_as_sf(coords = c('x', 'y')) %>%
                st_set_crs('EPSG:4326'), size = 0.05) +
      theme(legend.position = 'none')
  }
  
  n_distinct(d$poly_id)
  
  # using a subset of the data and polygons ----
  # need to fix the factor levels to ensure nbs and levels match
  poly_sub <- unique(d$poly_id)[1:1000] %>% factor(., levels = .)
  d_sub <- filter(d, poly_id %in% poly_sub)
  nbs_sub <- filter(ecoregions, poly_id %in% as.character(poly_sub)) %>%
    spdep::poly2nb(., row.names = .$poly_id, queen = TRUE)
  names(nbs_sub) <- attr(nbs_sub, 'region.id')
  d_sub$poly_id <- factor(d_sub$poly_id, levels = poly_sub)
  
  m <- bam(z ~ s(poly_id, bs = 'mrf', xt = list(nb = nbs_sub)) +
             s(eco_name, bs = 're'),
           family = gaussian(),
           data = d_sub,
           method = 'fREML',
           discrete = TRUE,
           drop.unused.levels = FALSE)
  summary(m)
  
  # plot the predictions for each shapefile
  ecoregions %>%
    filter(poly_id %in% as.character(poly_sub)) %>%
    rename(eco_name = ECO_NAME) %>%
    mutate(., mu_hat = predict(m, type = 'response', newdata = .)) %>%
    ggplot() +
    geom_sf(aes(fill = mu_hat)) +
    khroma::scale_fill_smoothrainbow()
}
