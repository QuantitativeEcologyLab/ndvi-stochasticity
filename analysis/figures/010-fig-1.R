library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('purrr')   # for functional programming
library('mgcv')    # for predicting from models
library('gratia')  # for working with GAMs
library('ggplot2') # for fancy plots
library('cowplot') # for fancy plots with multiple panels
source('analysis/figures/000-default-ggplot-theme.R')
source('analysis/figures/000-robinson-objects.R')

r_eco <- read_sf('data/ecoregions/ecoregions-polygons-edited.shp') %>%
  filter(group == 'Neotropic, Nearctic') %>%
  vect() %>%
  rasterize(rast('data/water-body-raster.tif'), field = 'biome')

shp_0 <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  st_geometry() %>%
  st_union() %>%
  st_as_sf()

shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  filter(group != 'Antarctic') %>%
  st_geometry() %>%
  st_union() %>%
  st_as_sf()
plot(shp, col = 'darkgreen')

r_0 <- rast('H:/GitHub/ndvi-stochasticity/data/avhrr-viirs-ndvi/raster-files/VIIRS-Land-v001/VIIRS-Land_v001_JP113C1_NOAA-20_20230101_c20240123205954.nc',
            lyr = 'NDVI')

if(FALSE) {
  d <-
    readr::read_csv('~/../Desktop/temp-preds.csv', show_col_types = FALSE) %>%
    rename(mu_hat_no_time = preds,
           mu_hat_2003_183 = preds_t) %>%
    mutate(denvar = map2_dbl(x, y, function(..x, ..y) {
      runif(1, 0,
            (cospi(..x / 180) * cospi(..y / 90) + 1.1) / 2.1 * 0.04)
    }))
} else {
  p_water <- rast('data/water-body-raster.tif')
  r_elev_m <- rast('data/elev-raster.tif')
  
  exclude <- smoooths(m_mu) %>%
    # c('s(prop_water)', 's(y,x)', 's(year)', 's(doy)', 's(elev_m)',
    #   'ti(year,doy)', 'ti(year,y,x)', 'ti(doy,y,x)', 'ti(year,doy,y,x)') %>%
    `[`(., which(grepl('doy', .) | grepl('year', .)))
  
  d <- r_0 %>%
    aggregate(4) %>%
    `values<-`(1) %>%
    mask(st_transform(shp, crs(.))) %>%
    as.data.frame(xy = TRUE) %>%
    as_tibble() %>%
    select(! NDVI) %>%
    mutate(prop_water = extract(p_water, tibble(x, y))[, 2],
           elev_m = extract(r_elev_m, tibble(x, y))[, 2]) %>%
    filter(prop_water < 0.4) %>%
    mutate(mu_hat = runif(n(), -1, 1),
           denvar = runif(n(), 0, 0.04))
  mutate(mu_hat = predict(m_mu, newdata = ., type = 'response',
                          se.fit = FALSE, exclude = exclude,
                          block.size = 1e6, discrete = FALSE),
         denvar = predict(m_s2, newdata = ., type = 'response',
                          se.fit = FALSE, exclude = exclude,
                          block.size = 1e6, discrete = FALSE))
  
  readr::write_csv(d, 'H:/temp/temp-preds.csv')
}

d_rob <- d %>%
  rast(crs = crs(r_0)) %>%
  project(robinson_crs) %>%
  as.data.frame(xy = TRUE) %>%
  as_tibble() %>%
  mutate(biome = extract(project(r_eco, robinson_crs), tibble(x, y))[, 2]) %>%
  filter(! (biome == 'Rock and Ice' & x > -6e6 & y > 0))

# Robinson projections of long-term mean NDVI and DENVar
fig_1 <- plot_grid(
  ggplot() +
    geom_sf(data = bounds, fill = 'white', color = 'black') +
    geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
    geom_raster(aes(x, y, fill = mu_hat_2003_183), d_rob) +
    geom_sf(data = shp_0, color = 'black', fill = 'transparent',
            lwd = 0.05) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_gradientn('Draft mean NDVI', colours = ndvi_pal,
                         limits = c(-1, 1)) +
    theme_void() +
    theme(legend.position = 'top', panel.grid = element_blank(),
          text = element_text(face = 'bold'), legend.key.width = rel(2)),
  ggplot() +
    geom_sf(data = bounds, fill = 'white', color = 'black') +
    geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
    geom_raster(aes(x, y, fill = denvar), d_rob) +
    geom_sf(data = shp_0, color = 'black', fill = 'transparent',
            lwd = 0.05) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_fill_davos(name = 'Fake DENVar', limits = c(0, NA), reverse = TRUE) +
    theme_void() +
    theme(legend.position = 'top', panel.grid = element_blank(),
          text = element_text(face = 'bold'), legend.key.width = rel(2)),
  labels = 'AUTO', ncol = 1)

ggsave('figures/fig-1.png', fig_1, width = 8.5, height = 10, dpi = 600,
       bg = 'white')

ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp_0, color = 'transparent', fill = 'grey') +
  geom_raster(aes(x, y, fill = biome == 'Rock and Ice' & ! is.na(biome)),
              d_rob) +
  geom_sf(data = shp_0, color = 'black', fill = 'transparent',
          lwd = 0.05) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual('Rock and Ice', values = c('grey30', 'cornflowerblue'),
                    labels = c('No', 'Yes')) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))
