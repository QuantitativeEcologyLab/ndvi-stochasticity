library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('purrr')   # for data wrangling
library('sf')      # for spatial data
library('terra')   # for rasters
library('ggplot2') # for fancy plots
library('khroma')  # for colorblind-friendly color palettes

# reference raster
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc',
            lyr = 'NDVI')

# get shapefile of each continent
ecoregions <- read_sf('data/wwf-ecoregions/wwf_terr_ecos.shp') %>%
  # remove 25 inland water polygons
  # great lakes are not among them, but there's no NDVI data for them
  filter(BIOME != 98) %>%
  # biomes cross-checked using QGIS
  mutate(
    # ecoregion < biome < realm < group
    ecoregion = ECO_NAME,
    biome = case_when(
      BIOME == 99 ~ 'Rock and Ice',
      BIOME == 11 ~ 'Tundra',
      BIOME == 6 ~ 'Boreal Forests/Taiga',
      BIOME == 10 ~ 'Montane Grasslands and Shrublands',
      BIOME == 8 ~ 'Temperate Grasslands, Savannas and Shrublands',
      BIOME == 4 ~ 'Temperate Broadleaf and Mixed Forests',
      BIOME == 5 ~ 'Temperate Coniferous Forests',
      BIOME == 12 ~ 'Mediterranean Forests, Woodlands and Scrub',
      BIOME == 14 ~ 'Mangroves',
      BIOME == 1 ~ 'Tropical and Subtropical Moist Broadleaf Forests',
      BIOME == 2 ~ 'Tropical and Subtropical Dry Broadleaf Forests',
      BIOME == 3 ~ 'Tropical and Subtropical Coniferous Forests',
      BIOME == 7 ~ 'Tropical and Subtropical Grasslands, Savannas and Shrublands',
      BIOME == 9 ~ 'Flooded Grasslands and Savannas',
      BIOME == 13 ~ 'Deserts and Xeric Shrublands') %>%
      factor(levels = c(
        'Rock and Ice', 'Tundra', 'Boreal Forests/Taiga',
        'Montane Grasslands and Shrublands',
        'Temperate Grasslands, Savannas and Shrublands',
        'Temperate Broadleaf and Mixed Forests',
        'Temperate Coniferous Forests',
        'Mediterranean Forests, Woodlands and Scrub',
        'Mangroves',
        'Tropical and Subtropical Moist Broadleaf Forests',
        'Tropical and Subtropical Dry Broadleaf Forests',
        'Tropical and Subtropical Coniferous Forests',
        'Tropical and Subtropical Grasslands, Savannas and Shrublands',
        'Flooded Grasslands and Savannas',
        'Deserts and Xeric Shrublands')),
    realm = case_when(REALM == 'AA' ~ 'Australasia',
                      REALM == 'AN' ~ 'Antarctic',
                      REALM == 'AT' ~ 'Afrotropic',
                      REALM == 'IM' ~ 'Indo-Malay',
                      REALM == 'NA' ~ 'Nearctic',
                      REALM == 'NT' ~ 'Neotropic',
                      REALM == 'OC' ~ 'Oceania',
                      REALM == 'PA' ~ 'Palearctic'),
    group = case_when(REALM %in% c('NA', 'NT') ~ 'americas',
                      REALM == 'AT' ~ 'afrotropic',
                      REALM == 'PA' ~ 'palearctic',
                      REALM %in% c('AA', 'IM', 'OC') ~ 'islands',
                      REALM == 'AN' ~ 'antarctica',
                      .default = 'unassigned')) %>%
  select(group, realm, biome, ecoregion, REALM) %>%
  st_make_valid() %>% # to drop duplicate vertices
  mutate(area_km2 = as.numeric(st_area(.)) / 1e6) %>%
  st_transform(crs = crs(r_0))

# add missing realms and groups by longitude
unassigned <- ecoregions %>%
  filter(group == 'unassigned') %>%
  bind_cols(.,
            st_coordinates(.) %>%
              data.frame() %>%
              group_by(L2) %>%
              summarize(min_x = min(X))) %>%
  select(! L2) %>%
  mutate(REALM = case_when(area_km2 > 9e6 ~ 'AN',
                           area_km2 > 1e6 ~ 'NA',
                           min_x > -30 ~ 'PA',
                           .default = 'NA'),
         realm_2 = case_when(REALM == 'AA' ~ 'Australasia',
                             REALM == 'AN' ~ 'Antarctic',
                             REALM == 'AT' ~ 'Afrotropic',
                             REALM == 'IM' ~ 'Indo-Malay',
                             REALM == 'NA' ~ 'Nearctic',
                             REALM == 'NT' ~ 'Neotropic',
                             REALM == 'OC' ~ 'Oceania',
                             REALM == 'PA' ~ 'Palearctic'),
         group_2 = case_when(REALM %in% c('NA', 'NT') ~ 'americas',
                             REALM == 'AT' ~ 'afrotropic',
                             REALM == 'PA' ~ 'palearctic',
                             REALM %in% c('AA', 'IM', 'OC') ~ 'islands',
                             REALM == 'AN' ~ 'antarctica',
                             .default = 'unassigned'))

filter(unassigned, group_2 == 'unassigned') # no polygons left unassigned

# make a shapefile of dissolved borders
shp_g <- ecoregions %>%
  summarize(geometry = st_union(geometry), .by = group) %>%
  mutate(area_km2 = st_area(.) / 1e6) %>%
  arrange(desc(area_km2)) %>%
  mutate(group = factor(group, levels = unique(group)))

plot(shp_g['group'])

# plot unassigned polygons by the new realm
ggplot() +
  geom_sf(data = shp_g, color = 'black', fill = 'grey') +
  geom_sf(aes(fill = realm_2), unassigned, color = NA) +
  scale_fill_highcontrast()

# ensure all polygons near the palearctic border are in the palearctic
ggplot() +
  geom_sf(data = shp_g, color = 'black', fill = 'grey') +
  geom_sf(data = unassigned, fill = 'red', color = NA) +
  geom_sf(data = shp_g, color = 'black', fill = NA) +
  coord_sf(xlim = c(65, 105), ylim = c(25, 45)) +
  scale_fill_highcontrast()

# label unassigned polygons
n_distinct(ecoregions$area_km2) == nrow(ecoregions) # can join by area

ecoregions <- left_join(ecoregions,
                        unassigned %>%
                          st_drop_geometry() %>%
                          select(group_2, realm_2, area_km2),
                        by = 'area_km2') %>%
  mutate(group = if_else(group == 'unassigned', group_2, group),
         realm = if_else(is.na(realm), realm_2, realm)) %>%
  select(! c(REALM, group_2, realm_2))

# ensure join is ok
ggplot() +
  geom_sf(aes(fill = group), ecoregions, color = 'black') +
  scale_fill_bright()

ggplot() +
  geom_sf(aes(fill = realm), ecoregions, color = 'black') +
  scale_fill_light()

# some realms have excessively small total areas
ecoregions %>%
  st_drop_geometry() %>%
  group_by(realm) %>%
  summarize(total = sum(area_km2)) %>%
  mutate(total = total / sum(total))

# groups have more reasonable total areas
ecoregions %>%
  st_drop_geometry() %>%
  group_by(group) %>%
  summarize(total = sum(area_km2)) %>%
  mutate(total = total / sum(total))

# update group names to include the different ecoregions
ecoregions <- mutate(ecoregions,
                     group = paste(unique(realm), collapse = ', '),
                     .by = group)
unique(ecoregions$group)

# update the shapefile of dissolved borders
shp_g <- ecoregions %>%
  summarize(geometry = st_union(geometry), .by = group) %>%
  mutate(area_km2 = st_area(.) / 1e6) %>%
  arrange(desc(area_km2)) %>%
  mutate(group = factor(group, levels = unique(group)))

# note: africa and arabian peninsula are split between groups 
ggplot() +
  geom_sf(aes(fill = group), shp_g, color = 'black') +
  scale_fill_bright() +
  theme_void() +
  theme(legend.position = 'top')

# save the final files ----
if(! dir.exists('data/ecoregions')) dir.create('data/ecoregions')
st_write(ecoregions, 'data/ecoregions/ecoregions-polygons.shp')

# edited in QGIS: made single african group, unified arabian peninsula
# made a shapefile of dissolved borders
ecoregions <- read_sf('data/ecoregions/ecoregions-polygons-edited.shp')
ecoregions

shp_g <- read_sf('data/ecoregions/groups-polygons.shp')
shp_g

# plot unassigned polygons by the new unified groups
ggplot() +
  geom_sf(aes(fill = group), data = shp_g, color = 'black') +
  scale_fill_bright() +
  theme(legend.position = 'top')
