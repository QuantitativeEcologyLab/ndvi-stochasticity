# lab PC: 128 GB of RAM, 12th Gen Intel(R) Core(TM) i7-12700K 3.60 GHz, 12 cores
# EME Linux: 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
#' *THE MODELS FIT IN THIS SCRIPT ARE FOR INITIAL TESTING ONLY*
#' *DO NOT USE THESE MODELS TO PREDICT FOR ANY PROJECTS*
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('purrr')   # for functional programming
library('furrr')   # for parallelized functional programming
library('mgcv')    # for Generalized Additive Models
library('gratia')  # for plotting GAMs
library('sf')      # for shapefiles
library('terra')   # for rasters
source('analysis/figures/000-default-ggplot-theme.R')

# shapefile of the world in Robinson projection
shp <- read_sf('data/ecoregions/groups-polygons.shp') %>%
  st_geometry() %>%
  st_union() %>%
  st_as_sf() %>%
  st_transform('ESRI:54030')

bounds <- tibble(x = c(-180:180, rep(180, 181), 180:-180, rep(-180, 181)),
                 y = c(rep(-90, 361), -90:90, rep(90, 361), 90:-90)) %>%
  st_as_sf(coords = c('x', 'y'), crs = 'EPSG:4326') %>%
  summarise(geometry = st_combine(geometry)) %>%
  st_cast('POLYGON') %>%
  st_transform(crs(shp))

# note: elevations are NA for abs(lat) > 85
# too memory-intensive to run on lab PC
if(file.exists('models/global-tests/models-with-100-days-only.rds')) {
  d <- readRDS('models/global-tests/models-with-100-days-only.rds')
} else {
  plan(multisession, workers = 4) # there are 4 different datasets
  
  d <-
    tibble(
      file = list.files('data/avhrr-viirs-ndvi/group-level-datasets',
                        pattern = '100-days', full.names = TRUE),
      group =
        gsub('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
             '', file) %>%
        gsub('-t-2-s-2-ndvi-data-100-days.rds', '', .),
      data = map(file, \(.fn) {
        readRDS(.fn) %>%
          na.omit() %>% #' missing slope and aspect at long = -180
          # add mean precipitation
          mutate(precip_m_day = extract(rast('data/precipitation-yearly-mean-m-day.tif'), tibble(x, y))[, 2])
      }, .progress = 'Importing data'),
      #' `k` values for `sos` smooth terms
      #' increasing `k` by >=100 gives the error:
      #'`Error in cbind(S, matrix(0, k - length(ind), length(ind))) :` 
      #'`  number of rows of matrices must match (see arg 2)`
      k_mu = c(2000, 900, 1300, 2000),
      k_s2 = c(2000, 900, 1300, 2000),
      m_mu = future_map2(data, k_mu, function(.d, .k) {
        bam(ndvi_aggr ~
              s(prop_water, bs = 'cr', k = 5) + # water biases to < 0
              s(sqrt(elev_m), bs = 'cr', k = 5) + # elevation effect
              s(y, x, bs = 'sos', k = .k) + # smooth of space
              s(aspect_deg, bs = 'cc', k = 10) + # aspect is circular
              s(sqrt(sqrt(slope_deg)), k = 6) + # very few values > 4
              s(log(precip_m_day), bs = 'cr', k = 6) +
              ti(sqrt(elev_m), y, x, bs = c('cr', 'sos'), d = c(1, 2),
                 k = c(6, 100)),
            family = gaussian(),
            data = .d, method = 'fREML', discrete = TRUE, nthreads = 12)
      }, .progress = TRUE, .options = furrr_options(seed = NULL)),
      sos_edf_mu = map_dbl(m_mu, function(.m) {
        sum(.m$edf[which(grepl('s\\(y,x\\)', names(.m$edf)))])
      })) %>%
    #' can't edit `data` column in `tibble()`
    # add predicted values and squared residuals
    mutate(
      data = future_map2(data, m_mu, function(.d, .m) {
        mutate(.d,
               mu_hat = predict(.m, discrete = FALSE, block.size = 5e5,
                                newdata.guaranteed = TRUE),
               e = ndvi_aggr - mu_hat,
               e2 = e^2)
      }, .progress = TRUE, .options = furrr_options(seed = NULL)),
      m_s2 = future_map2(data, k_s2, function(.d, .k) {
        # s2_hat = E((Y - mu_hat)^2) = E(e^2)
        bam(e2 ~
              s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
              s(sqrt(elev_m), bs = 'cr', k = 5) + # elevation effect
              s(y, x, bs = 'sos', k = .k) + # smooth of space
              s(aspect_deg, bs = 'cc', k = 10) + # aspect is circular
              s(sqrt(sqrt(slope_deg)), k = 6) + # very few values > 4
              s(log(precip_m_day), bs = 'cr', k = 6) +
              ti(sqrt(elev_m), y, x, bs = c('cr', 'sos'), d = c(1, 2),
                 k = c(6, 100)) +
              s(mu_hat, bs = 'cr', k = 10), # mu and s2 are correlated
            family = gaussian(),
            data = .d, method = 'fREML', discrete = TRUE, nthreads = 12)
      }, .progress = TRUE, .options = furrr_options(seed = NULL)),
      sos_edf_s2 = map_dbl(m_s2, function(.m) {
        sum(.m$edf[which(grepl('s\\(y,x\\)', names(.m$edf)))])
      }),
      group_label = stringr::str_to_title(group) %>%
        gsub('-', ', ', .) %>%
        gsub('Indo, Malay', 'Indo-Malay', .))
  
  plan(sequential)
  
  saveRDS(d, 'models/global-tests/models-with-100-days-only.rds')
}

# clear the text file, if it exists
sink('models/global-tests/100-days-model-summaries.txt', append = FALSE)
cat('')
sink()

for(i in 1:nrow(d)) {
  cat('Working on model', i, 'of', nrow(d), ' (', d$group_label[i], ')', '\n')
#' `draw()` does not plot much below 65 degrees S
#' *note:* x and y axes are flipped for the `sos` smooth
png(paste0('figures/global-tests/', d$group[i], '-test-models.png'),
    width = 6, height = 24, units = 'in', res = 300, bg = 'white')
layout(matrix(1:16, ncol = 2, byrow = FALSE))
plot(d$m_mu[[i]], rug = FALSE, scheme = c(1, 1, 5, 1, 1, 1, 3),
     n2 = 75, too.far = 0.05, scale = 0)
plot.new() # to fill the bottom plot of the first column
plot(d$m_s2[[i]], rug = FALSE, scheme = c(1, 1, 5, 1, 1, 1, 3, 1),
     n2 = 75, too.far = 0.05, scale = 0)
dev.off()

sink('models/global-tests/100-days-model-summaries.txt', append = TRUE)
cat('region: ', d$group[i], '; parameter: mean\n', sep = '')
print(summary(d$m_mu[[i]]))
cat('\n- - - - - - - - - - - - - -\n', sep = '')
print(summary(d$m_s2[[i]]))
cat('\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n', sep = '')
sink() # close connection to text file
}; rm(i)

# check pixel-level mean residuals
if(file.exists('models/global-tests/data-summary-from-models-with-100-days-only.rds')) {
  d_summ <- readRDS('models/global-tests/data-summary-from-models-with-100-days-only.rds')
} else {
  d_summ <- d %>%
    select(group, group_label, data) %>%
    unnest(data) %>%
    select(group, group_label, e, mu_hat, x, y, prop_water, elev_m,
           slope_deg, aspect_deg, precip_m_day) %>%
    summarize(n = n(),
              mu_hat = mean(mu_hat, na.rm = TRUE),
              mean_e = mean(e, na.rm = TRUE),
              prop_water = mean(prop_water, na.rm = TRUE),
              elev_m = mean(elev_m, na.rm = TRUE),
              slope_deg = mean(slope_deg, na.rm = TRUE),
              aspect_deg = mean(aspect_deg, na.rm = TRUE),
              precip_m_day = mean(precip_m_day, na.rm = TRUE),
              .by = c(group, group_label, x, y)) %>%
    nest(dat = ! c(group, group_label)) %>%
    mutate(dat = map2(dat, group, \(.d, .g) {
      .m <- d$m_s2[[which(d$group == .g)]]
      mutate(.d,
             s2_hat = predict.bam(.m, newdata = .d, type = 'response',
                                  se.fit = FALSE, newdata.guaranteed = TRUE,
                                  discrete = FALSE)) %>%
        return()
    })) %>%
    unnest(dat)
  saveRDS(d_summ, 'models/global-tests/data-summary-from-models-with-100-days-only.rds')
  
  if(FALSE) {
    hist(d_summ$slope_deg)
    mean(d_summ$slope_deg > 4, na.rm = TRUE)
    
    layout(1:2)
    hist(d_summ$slope_deg^(1/4)) # much less kurtosis after transformation
    hist(d_summ$aspect_deg) # very close to uniform distribution
    layout(1)
  }
}

# import the first raster
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/AVHRR-Land_v006/N07_AVH13C1/N07_AVH13C1.A1981175.006.2022270161458.nc',
            lyr = 'NDVI')

#' need to split by group to ensure raster is not too large
d_summ_rob <-
  d_summ %>%
  select(x, y, n, mu_hat, mean_e, s2_hat, group_label) %>%
  nest(g_data = ! c(group_label)) %>%
  mutate(g_data = map(g_data, \(.d) {
    .d %>%
      rast(type = 'xyz', crs = crs(r_0)) %>%
      project(crs(shp)) %>%
      as.data.frame(xy = TRUE)
  })) %>%
  unnest(g_data) %>%
  mutate(mean_e = case_when(mean_e > 0.2 ~ 0.2,
                            mean_e < -0.2 ~ -0.2,
                            .default = mean_e),
         s2_hat = if_else(s2_hat < 0, 4.406384e-09, s2_hat))

#' single `geom_raster` layer gives: `cannot allocate vector of size 2438889.4 Gb`
p_map_n <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp, fill = 'grey30') +
  geom_raster(aes(x, y, fill = n),
              filter(d_summ_rob, group_label == 'Africa')) +
  geom_raster(aes(x, y, fill = n),
              filter(d_summ_rob, group_label == 'Indo-Malay, Oceania, Australasia')) +
  geom_raster(aes(x, y, fill = n),
              filter(d_summ_rob, group_label == 'Neotropic, Nearctic')) +
  geom_raster(aes(x, y, fill = n),
              filter(d_summ_rob, group_label == 'Palearctic')) +
  scale_fill_bamako(name = 'Number of non-NA cells', reverse = TRUE) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))

p_map_mu <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp, fill = 'grey30') +
  geom_raster(aes(x, y, fill = mu_hat),
              filter(d_summ_rob, group_label == 'Africa')) +
  geom_raster(aes(x, y, fill = mu_hat),
              filter(d_summ_rob, group_label == 'Indo-Malay, Oceania, Australasia')) +
  geom_raster(aes(x, y, fill = mu_hat),
              filter(d_summ_rob, group_label == 'Neotropic, Nearctic')) +
  geom_raster(aes(x, y, fill = mu_hat),
              filter(d_summ_rob, group_label == 'Palearctic')) +
  scale_fill_gradientn('Mean NDVI', colors = ndvi_pal,
                       limits = c(-1, 1)) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))

p_map_mean_e <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp, fill = 'grey30') +
  geom_raster(aes(x, y, fill = mean_e),
              filter(d_summ_rob, group_label == 'Africa')) +
  geom_raster(aes(x, y, fill = mean_e),
              filter(d_summ_rob, group_label == 'Indo-Malay, Oceania, Australasia')) +
  geom_raster(aes(x, y, fill = mean_e),
              filter(d_summ_rob, group_label == 'Neotropic, Nearctic')) +
  geom_raster(aes(x, y, fill = mean_e),
              filter(d_summ_rob, group_label == 'Palearctic')) +
  scale_fill_distiller('Mean residual, capped', type = 'div', palette = 5,
                       limits = c(-0.2, 0.2)) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))

p_map_s2 <-
  ggplot() +
  geom_sf(data = bounds, fill = 'white', color = 'black') +
  geom_sf(data = shp, fill = 'grey30') +
  geom_raster(aes(x, y, fill = s2_hat),
              filter(d_summ_rob, group_label == 'Africa')) +
  geom_raster(aes(x, y, fill = s2_hat),
              filter(d_summ_rob, group_label == 'Indo-Malay, Oceania, Australasia')) +
  geom_raster(aes(x, y, fill = s2_hat),
              filter(d_summ_rob, group_label == 'Neotropic, Nearctic')) +
  geom_raster(aes(x, y, fill = s2_hat),
              filter(d_summ_rob, group_label == 'Palearctic')) +
  scale_fill_davos(name = 'Variance in NDVI', limits = c(0, 0.1)) +
  theme_void() +
  theme(legend.position = 'top', panel.grid = element_blank(),
        text = element_text(face = 'bold'), legend.key.width = rel(2))

# proportion of estimated variances not in the scale's range
mean(d_summ_rob$s2_hat == 4.406384e-09) # 4.6%
mean(d_summ_rob$s2_hat > 0.1) # 0.01%

p_map <- cowplot::plot_grid(p_map_n, p_map_mu, p_map_mean_e, p_map_s2,
                            labels = 'AUTO', ncol = 2)

ggsave('figures/global-tests/100-days-summary-map.png', p_map,
       width = 20, height = 12, units = 'in', dpi = 300, bg = 'white')
