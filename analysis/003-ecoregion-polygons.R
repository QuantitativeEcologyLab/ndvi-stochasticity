library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('purrr')   # for data wrangling
library('sf')      # for spatial data
library('terra')   # for rasters
library('spdep')   # for finding neighbors
library('ggplot2') # for fancy plots
sf_use_s2(FALSE)

# get shapefile of each continent
ecoregions <- st_read('data/world-ecosystems/data/commondata/data0/tnc_terr_ecoregions.shp') %>%
  mutate(group = case_when(WWF_REALM %in% c('NA', 'NT') ~ 'americas',
                           WWF_REALM == 'AT' ~ 'afrotropic',
                           WWF_REALM == 'PA' ~ 'paleoartctic',
                           WWF_REALM %in% c('AA', 'IM', 'OC') ~ 'islands',
                           WWF_REALM == 'AN' ~ 'antarctica',
                           .default = 'unassigned')) %>%
  select(group, ECO_NAME, ECO_NUM, WWF_REALM, WWF_MHTNAM) %>%
  st_make_valid() %>% # to drop duplicate vertices
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6) %>%
  filter(WWF_MHTNAM != 'Inland Water') %>% # rm 17 inland water polygons
  st_transform(crs = crs(rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v005_AVH13C1_NOAA-07_19810624_c20170610041337.nc',
                              lyr = 'NDVI')))

if(FALSE) { # some exploratory plots and checks
  # each ecoregion has a unique name, but many ecoregions are multipolygons
  unique(ecoregions$ECO_NAME)[1:10]
  plot(st_geometry(ecoregions), col = factor(ecoregions$ECO_NAME),
       axes = TRUE)
  
  # by group ----
  # find total area by group
  ecoregions %>%
    st_drop_geometry() %>%
    group_by(group) %>%
    summarize(total = sum(area_km2)) %>%
    mutate(total = total / sum(total))
  
  ggplot(ecoregions) +
    facet_wrap(~ group) +
    geom_sf(fill = 'red3', color = 'black') +
    theme_classic()
  
  # by WWF realm ----
  # find total area by WWF_REALM
  ecoregions %>%
    st_drop_geometry() %>%
    group_by(WWF_REALM) %>%
    summarize(total = sum(area_km2)) %>%
    mutate(total = total / sum(total))
  
  ggplot(ecoregions) +
    facet_wrap(~ WWF_REALM) +
    geom_sf(fill = 'red3', color = 'black') +
    theme_classic()
}

# save the final files ----
if(! dir.exists('data/ecoregions')) dir.create('data/ecoregions')
st_write(ecoregions, 'data/ecoregions/ecoregions-polygons.shp')
