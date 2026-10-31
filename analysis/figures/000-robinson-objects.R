library('sf')

robinson_crs <- 'ESRI:54030'

bounds <-
  dplyr::tibble(x = c(-180:180, rep(180, 181), 180:-180, rep(-180, 181)),
                y = c(rep(-90, 361), -90:90, rep(90, 361), 90:-90)) %>%
  st_as_sf(coords = c('x', 'y'), crs = 'EPSG:4326') %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast('POLYGON') %>%
  st_transform(robinson_crs)

shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  st_geometry() %>%
  st_as_sf() %>%
  st_union() %>%
  st_as_sf()

shp_rob <- st_transform(shp, robinson_crs)
