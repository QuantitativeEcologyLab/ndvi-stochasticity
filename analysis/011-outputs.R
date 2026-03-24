library('sf')     # for shapefiles
library('mgcv')   # for predicting from models
library('gratia') # for selecting model terms
library('terra')  # for working with rasters
library('dplyr')  # for data wrangling
library('purrr')  # for functional programming
library('tidyr')  # for data wrangling

GROUPS <- c('neotropic-nearctic', 'africa',
            'indo-malay-oceania-australasia', 'palearctic')

for(GROUP in GROUPS) {
  shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
    mutate(group = tolower(gsub(' ', '-', group))) %>%
    filter(group == GROUP) %>%
    select(group)
  
  r_prop_water <- rast('data/water-body-raster.tif') %>%
    ifel(. > 0.4, NA, .) %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE) %>%
    aggregate(4, na.rm = TRUE)
  
  
  r_elev <- rast('data/elev-raster.tif') %>%
    crop(st_transform(shp, crs = crs(.)), mask = TRUE) %>%
    project(shp)
  
  layout(1:2)
  plot(r_prop_water)
  plot(r_elev)
  layout(1)
  
  m_mu <- qread(list.files(path = 'models/global-models',
                           pattern = paste0('bam-mean-ndvi-', GROUP),
                           full.names = TRUE),
                nthreads = NCORES)
  m_s2 <- qread(list.files(path = 'models/global-models',
                           pattern = paste0('bam-variance-ndvi-', GROUP),
                           full.names = TRUE),
                nthreads = NCORES)
  
  #' exclude all terms with `year` or `doy`:
  #' `(Intercept)`: included
  #' `s(prop_water)`: included
  #' `s(y,x)`: included
  #' `s(year)`: excluded
  #' `s(doy)`: excluded
  #' `s(elev_m)`: included
  #' `s(mu_hat)`: included, not present in mean model
  #' `ti(year,doy)`: excluded
  #' `ti(year,y,x)`: excluded
  #' `ti(doy,y,x)`: excluded
  #' `ti(year,doy,y,x)`: excluded
  EXCLUDE <- smooths(m_s2)[grepl('year', smooths(m_s2)) |
                             grepl('doy', smooths(m_s2))]
  
  # runs within ~10 minutes per group
  preds <- 
    list(r_prop_water, project(r_elev, r_prop_water)) %>%
    rast() %>%
    as.data.frame(xy = TRUE, na.rm = TRUE) %>%
    as_tibble() %>%
    rename(prop_water = layer,
           elev_m = file37905ccd4e9d) %>%
    mutate(year = 0, doy = 0) %>% # excluded, so values don't matter
    mutate(mu_hat = predict(m_mu, newdata = ., discrete = FALSE,
                            type = 'response', se.fit = FALSE,
                            exclude = EXCLUDE)) %>%
    # estimated mean is necessary for the estimated variance
    mutate(s2_hat = predict.bam(m_s2, newdata = ., discrete = FALSE,
                                type = 'response', se.fit = FALSE,
                                exclude = EXCLUDE))
  preds
  
  saveRDS(preds, paste0('output/long-term-preds-', GROUP, '.rds'))
  
  rm(preds, m_mu, m_s2, r_prop_water, r_elev, shp)
  gc()
}

# EME linux has issues with saving the rasters
# (probably a GDAL installation issue after the last update)
extent <- read_sf('data/ecoregions/ecoregions-polygons-edited.shp') %>%
  filter(group != 'Antarctic') %>%
  filter(biome != 'Rock and Ice') %>%
  # predictions for Kalaallit Nunaat are poor because of data scarcity
  filter(! grepl('Kalaallit Nunaat', ecoregion))

rasters <-
  tibble(file_name = list.files(path = 'H:/GitHub/ndvi-stochasticity/output',
                                pattern = 'long-term-preds-.*.rds',
                                full.names = TRUE),
         group = gsub('.*-preds-', '', file_name)) %>%
  mutate(group = gsub('.rds', '', group),
         est = map2(file_name, group, function(.fn, .g) {
           readRDS(.fn) %>%
             select(x, y, mu_hat, s2_hat) %>%
             mutate(mu_hat = if_else(mu_hat < -1, -1, mu_hat),
                    s2_hat = if_else(s2_hat < 0, 0, s2_hat)) %>%
             rast() %>%
             `crs<-`('EPSG:4326') %>%
             crop(st_transform(extent, st_crs(.))) %>%
             project(r_prop_water)
         })) %>%
  pull(est) %>%
  sds() %>%
  app(fun = 'sum', na.rm = TRUE) %>%
  mask(st_transform(extent, st_crs(.))) %>%
  ifel(init(., 'y') < 70, ., NA) %>%
  `names<-`(c('mu_hat', 's2_hat'))

par(mfrow = c(2, 1))
plot(rasters)

writeRaster(rasters, paste0('output/long-term-preds.tif'))
