# lab PC: 12th Gen Intel(R) Core(TM) i7-12700K 3.60 GHz, 128 GB of RAM
#' *THE MODELS FIT IN THIS SCRIPT ARE FOR INITIAL TESTING ONLY*
#' *DO NOT USE THESE MODELS TO PREDICT*
library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('purrr')   # for functional programming
library('furrr')   # for parallelized functional programming
library('mgcv')    # for Generalized Additive Models
library('gratia')  # for plotting GAMs
library('sf')      # for shapefiles
source('analysis/figures/000-default-ggplot-theme.R')

# note: elevations are NA for abs(lat) > 85
if(file.exists('models/global-tests/models-with-first-100-days-only.rds')) {
  d <- readRDS('models/global-tests/models-with-first-100-days-only.rds')
} else {
  d <-
    tibble(file = list.files('data/avhrr-viirs-ndvi/group-level-datasets',
                             pattern = 'first-100-days', full.names = TRUE),
           group =
             gsub('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-',
                  '', file) %>%
             gsub('-t-2-s-2-ndvi-data-first-100-days.rds', '', .),
           # dropping NA elevations in Antarctica
           data = map(file, \(fn) na.omit(readRDS(fn)),
                      .progress = 'Importing data'),
           k_mu = c(500, 500, 500, 500, 500),
           k_s2 = c(500, 500, 500, 500, 500),
           m_mu = map2(data, k_mu, function(.d, .k) {
             bam(ndvi_aggr ~
                   s(prop_water, bs = 'cr', k = 5) + # water biases to < 0
                   s(elev_m, bs = 'cr', k = 5) + # elevation effect
                   s(y, x, bs = 'sos', k = .k),# + # smooth of space
                 family = gaussian(), # similar and faster results to betar()
                 data = .d, method = 'fREML', discrete = TRUE)
           }, .progress = 'Fitting mean models'),
           sos_edf_mu = map_dbl(m_mu, function(.m) {
             sum(.m$edf[which(grepl('s\\(y,x\\)', names(.m$edf)))])
           })) %>%
    mutate(
      # add predicted values and squared residuals
      data = map2(data, m_mu, function(.d, .m) {
        mutate(.d,
               mu_hat = predict(.m, discrete = TRUE), # to keep things fast
               e = (ndvi_aggr - mu_hat)^2,
               e2 = e^2)
      }, .progress = 'Adding residuals'),
      m_s2 = map2(data, k_s2, function(.d, .k) {
        bam(e2 ~ # s2_hat = E((Y - mu_hat)^2) = E(e^2)
              s(prop_water, bs = 'cr', k = 5) + # water biases towards < 0
              s(elev_m, bs = 'cr', k = 5) + # elevation effect
              s(y, x, bs = 'sos', k = .k), # smooth of space
            family = gaussian(),
            data = .d, method = 'fREML', discrete = TRUE)
      }, .progress = 'Fitting variance models'),
      sos_edf_s2 = map_dbl(m_s2, function(.m) {
        sum(.m$edf[which(grepl('s\\(y,x\\)', names(.m$edf)))])
      }))
  
  saveRDS(d, 'models/global-tests/models-with-first-100-days-only.rds')
}

for(i in 1:nrow(d)) {
  cat('Working on model', i, 'of', nrow(d), '\n')
  #' `draw()` does not plot much below 65 degrees S
  #' *note:* x and y axes are flipped for the `sos` smooth
  png(paste0('figures/global-tests/', d$group[i], '-test-model-mu.png'),
      width = 12, height = 12, units = 'in', res = 300, bg = 'white')
  layout(matrix(1:6, ncol = 2, byrow = FALSE))
  plot(d$m_mu[[i]], rug = FALSE, scheme = c(1, 1, 5), n2 = 75,
       too.far = 0.05, scale = 0)
  plot(d$m_s2[[i]], rug = FALSE, scheme = c(1, 1, 5), n2 = 75,
       too.far = 0.05, scale = 0)
  dev.off()
  
  sink('models/global-tests/first-100-days-model-summaries.txt', append = TRUE)
  cat('region: ', d$group[i], '; parameter: mean\n', sep = '')
  print(summary(d$m_mu[[i]]))
  cat('\n--------------------------\n', sep = '')
  print(summary(d$m_s2[[i]]))
  cat('\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =\n\n', sep = '')
  sink() # close connection to text file
}; rm(i)

# check pixel-level mean residuals
mean_e <- d %>%
  select(group, data) %>%
  unnest(data) %>%
  select(group, x, y, e) %>%
  summarize(mean_e = mean(e, na.rm = TRUE), .by = c(group, x, y))

p_hist <-
  ggplot(mean_e) +
  facet_wrap(~ group, scales = 'free', ncol = 2) +
  geom_histogram(aes(mean_e), binwidth = 0.01, center = 0.005,
                 fill = 'grey', color = 'black') +
  xlab('Mean pixel-level residual') +
  scale_y_log10(expression(bold(Count~(log['10']~scale))))

ggsave('figures/global-tests/first-100-days-mean-residuals-histogram.png',
       p_hist, width = 8, height = 8, units = 'in', dpi = 300, bg = 'white')

shp <- read_sf('data/ecoregions/groups-polygons.shp')

p_map <-
  ggplot(mean_e) +
  geom_sf(data = shp, fill = 'grey30') +
  geom_raster(aes(x, y, fill = mean_e)) +
  scale_fill_distiller('Mean residual', limits = c(-0.1, 0.1),
                       type = 'div', palette = 5) +
  labs(x = NULL, y = NULL)

ggsave('figures/global-tests/first-100-days-mean-residuals-map.png',
       p_map, width = 12, height = 6, units = 'in', dpi = 300, bg = 'white')
