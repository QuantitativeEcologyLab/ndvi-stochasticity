# specs of the EME Linux:
# 2.2 TB RAM, Intel Xeon Platinum 8462Y+ processor, 64 cores
#' some code wrappted in `if(FALSE){...}` because it takes too long to run
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
source('functions/scale-ndvi.R')
source('functions/ndvi-palette.R')
source('functions/decode_qa.R')    # for cleaning rasters
source('functions/bit_to_int.R')   # for cleaning rasters
source('functions/is_flagged.R') # for cleaning rasters
source('analysis/figures/000-default-ggplot-theme.R')
source('functions/get_legend.R') # get_legend() from cowplot v. 1.1.3 fails
source('functions/nbs_from_rast.R') # gives a list of neighboring cells

ecoregions <- read_sf('data/ecoregions/ecoregions-polygons.shp')
taiga <- filter(ecoregions, ecoregion == 'Northwest Territories taiga')

# some very large NDVI values at the edge of the preliminary cleaning ----
plot(st_geometry(ecoregions))
plot(st_geometry(taiga), add = TRUE, col = 'red')

list.files('data/avhrr-viirs-ndvi/raster-files',
           'AVHRR-Land_v005_AVH13C1_NOAA-07_1982011',
           full.names = TRUE) %>%
  map(\(x) rast(x, lyr = 'NDVI')) %>%
  rast() %>%
  crop(taiga, mask = TRUE) %>%
  plot()

# import raster for 2024-01-10
r_0 <- rast('data/avhrr-viirs-ndvi/raster-files/VIIRS-Land_v001_JP113C1_NOAA-20_20240110_c20240124141257.nc') %>%
  crop(ecoregions, mask = TRUE)
plot(r_0)

decoded_0 <- decode_qa(as.numeric(values(r_0$QA)))

# make rasters for each QA parameter
for(cn in colnames(decoded_0)[3:ncol(decoded_0)]) {
  r_0[[cn]] <- r_0$QA
  values(r_0[[cn]]) <- decoded_0[[cn]]
  r_0[[cn]] <- mask(r_0[[cn]], r_0$QA)
}
plot(r_0)

d_0 <- tibble(as.data.frame(r_0, na.rm = TRUE, xy = TRUE))
range(d_0$y)
summary(filter(d_0, NDVI > 0.2)) # max lat is too high for such high NDVI
plot(NDVI ~ y, filter(d_0, y > 60))

# cloud and snow/ice flags aren't very reliable
plot(r_0$cloud_flag)
d_0 %>%
  ggplot() +
  facet_grid(snow_ice ~ land_type) +
  geom_boxplot(aes(cloud_state, NDVI), alpha = 0.01)

d_0 %>%
  summarise(n = n(),
            median = median(NDVI, na.rm = TRUE),
            q_0.75 = quantile(NDVI, 0.75, na.rm = TRUE),
            upr_IQR = quantile(NDVI, 0.75) + 1.5 * IQR(NDVI, na.rm = TRUE),
            .by = c(snow_ice, cloud_state))
rm(r_0)

# first year of data, every 10 days
r_0 <- map(list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                      pattern = '.nc',
                      full.names = TRUE)[seq(1, 365, by = 10)],
           \(fn) rast(fn, lyr = 'NDVI')) %>%
  rast()
plot(r_0)

# make some rasters of summary stats every 10 rasters for each month ----
# max(NDVI); runs in ~ 10 minutes on EME Linux ----
if(! file.exists('data/avhrr-viirs-ndvi/max-raster-jan.tif')) {
  for(i in 1:12) {
    print(i)
    month <- tolower(month.abb[i])
    
    r_max <-
      map(list.files(
        path = 'data/avhrr-viirs-ndvi/raster-files',
        pattern = if_else(i < 10, paste0('0', i), as.character(i)) %>%
          paste0(, '.._.*\\.nc'),
        full.names = TRUE) %>%
          `[`(., seq(1, length(.), by = 10)), # take one every 10 days
        \(.fn) .ndvi <- rast(.fn, lyr = 'NDVI')) %>%
      rast() %>%
      max(na.rm = TRUE)
    
    writeRaster(r_max, paste0('data/avhrr-viirs-ndvi/max-raster-',
                              month, '.tif'), overwrite = TRUE)
    rm(i, month, r_max)
  }
}

# plot rasters of max values
png('figures/global-diagnostics/monthly-max-values.png',
    width = 15, height = 5.5, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12) {
  plot(rast(paste0('data/avhrr-viirs-ndvi/max-raster-',
                   tolower(month.abb[i]),
                   '.tif')),
       main = month.name[i])
  rm(i)
}
dev.off()

# histograms of max values at high latitudes
png('figures/global-diagnostics/monthly-max-values-hist.png',
    width = 15, height = 15, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12) {
  rast(paste0('data/avhrr-viirs-ndvi/max-raster-', tolower(month.abb)[i],
              '.tif')) %>%
    as.data.frame(xy = TRUE) %>%
    filter(max > 0.2, y > 50) %>%
    pull(y) %>%
    hist(main = paste('Cells with max(NDVI) > 0.2 in', month.name[i]),
         breaks = 1000, xlim = c(50, 85), ylim = c(0, 6e3),
         xlab = 'Latitude')
  abline(v = 58, col = 'red')
  rm(i)
}
dev.off()

# dropping top 2 degrees ----
# we could clean the rasters by simply dropping some top rows in the
# problematic months, but this also drops good data, too, and data is
# sparse in winter at high latitudes. Two degrees also does not seem
# to be enough (see February in the figure below).
# I (Stefano) couldn't find any QA flags associated with the high-lat bands
if(! file.exists('data/avhrr-viirs-ndvi/max-raster-dec-without-top-2-degrees.tif')) {
  for(i in 1:12) {
    print(i)
    month <- tolower(month.abb[i])
    month_char <- as.character(i)
    month_char <- if_else(nchar(month_char) == 2, month_char,
                          paste0('0', month_char))
    
    r_max_clean <-
      map(list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                     pattern = paste0(month_char, '.._.*\\.nc'),
                     full.names = TRUE) %>%
            `[`(., seq(1, length(.), by = 10)),
          \(.fn) {
            .ndvi <- rast(.fn, lyr = 'NDVI')
            
            # drop first few rows of data since they often include errors
            #' rows 50-90 N: `nrow(z) / (90*2) * (90-50) = 800`
            row_1 <- floor(which(is.finite(.ndvi[1:800, ]$NDVI))[1] /
                             ncol(.ndvi))
            .ndvi[row_1 + 0:39, ] <- NA
            return(.ndvi)
            
            # THERE IS NO QA FLAG THAT REMOVES THE HIGH-NDVI BANDS
            # .qa <- rast(.fn, lyr = 'QA')
            # 
            # to drop values based on QA flags
            # ifel(
            # # NOTE: bits are written left to right in R
            # # drop "confidently cloudy" pixels
            # (is_flagged(.qa, flag_position = 0) &
            #    is_flagged(.qa, flag_position = 1)) |
            #   # drop pixels with "cloud shadow"
            #   is_flagged(.qa, flag_position = 2) |
            #   # flags 3-5 (land/water) are dropped using spatial masking
            #   # # drop pixels with aerosol quality that is "poor" ("! OK")
            #   # ! is_flagged(.qa, flag_position = 6) |
            #   # flags 7, 11-14 are unused
            #   # drop pixels with thin cirrus reflective clouds
            #   is_flagged(.qa, flag_position = 8) |
            #   # drop pixels with thin cirrus emissive clouds
            #   is_flagged(.qa, flag_position = 9) |
            #   # drop pixels with cloud flag
            #   is_flagged(.qa, flag_position = 10)
            # # snow/ice flag is ok as both 0 or 1
            # ,
            # yes = NA, no = .ndvi) %>%
            # return()
          }) %>%
      rast() %>%
      max(na.rm = TRUE)
    
    writeRaster(
      x = r_max_clean,
      filename = paste0('data/avhrr-viirs-ndvi/max-raster-',
                        month, '-without-top-2-degrees.tif'),
      overwrite = TRUE)
  }
}

png('figures/global-diagnostics/monthly-max-values-without-top-2-degrees.png',
    width = 15, height = 5.5, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12){
  plot(rast(paste0('data/avhrr-viirs-ndvi/max-raster-', tolower(month.abb)[i],
                   '-without-top-2-degrees.tif')),
       main = month.name[i])
  rm(i)
}
dev.off()

# histograms of max values at high latitudes
png('figures/global-diagnostics/monthly-max-values-without-top-2-degrees-hist.png',
    width = 15, height = 15, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12) {
  rast(paste0('data/avhrr-viirs-ndvi/max-raster-', tolower(month.abb)[i],
              '-without-top-2-degrees.tif')) %>%
    as.data.frame(xy = TRUE) %>%
    filter(max > 0.2, y > 50) %>%
    pull(y) %>%
    hist(main = paste('Cells with max(NDVI) > 0.2 in', month.name[i]),
         breaks = 1000, xlim = c(50, 85), ylim = c(0, 6e3),
         xlab = 'Latitude')
  abline(v = 58, col = 'red')
  rm(i)
}
dev.off()

# filtering for ndvi < 90th percentile doesn't fix the issue ----
# the resulting rasters still contain problematic values
if(file.exists('data/avhrr-viirs-ndvi/90-percentile-raster-jan.tif')) {
  r_jan_90 <- rast('data/avhrr-viirs-ndvi/90-percentile-raster-jan.tif')
  plot(r_jan_90)
} else {
  r_jan_90 <-
    map(list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                   pattern = '01.._.*\\.nc', full.names = TRUE) %>%
          `[`(., seq(1, length(.), by = 10)),
        \(.fn) rast(.fn, lyr = 'NDVI')) %>%
    rast() %>%
    quantile(0.90, na.rm = TRUE)
  plot(r_jan_90)
  writeRaster(r_jan_90, 'data/avhrr-viirs-ndvi/90-percentile-raster-jan.tif')
}

# no clear breaks in the values
r_jan_90_top <-
  crop(r_jan_90,
       vect(tibble(x = c(-179, 180), y = c(59, 70)), geom = c('x','y')) %>%
         st_bbox() %>%
         st_as_sfc() %>%
         st_as_sf() %>%
         st_set_crs('EPSG:4326') %>%
         st_transform(crs(r_jan_90)), mask = TRUE)
plot(r_jan_90_top)
hist(r_jan_90_top, breaks = 1000, ylim = c(0, 1000))
# most values are 0
sort(table(round(values(r_jan_90_top), 3)), decreasing = TRUE)[1:3]

# july doesn't have the same issues (ignoring antarctica) ----
# despite it being winter in the southern hemisphere
if(file.exists('data/avhrr-viirs-ndvi/max-raster-jul.tif')) {
  plot(rast('data/avhrr-viirs-ndvi/max-raster-jul.tif'))
} else {
  r_jul_max <-
    map(list.files(path = 'data/avhrr-viirs-ndvi/raster-files',
                   pattern = '07.._.*\\.nc', full.names = TRUE) %>%
          `[`(., seq(1, length(.), by = 10)),
        \(.fn) rast(.fn, lyr = 'NDVI')) %>%
    rast() %>%
    max(na.rm = TRUE)
  plot(r_jul_max)
  writeRaster(r_jul_max, 'data/avhrr-viirs-ndvi/max-raster-jul.tif')
}

# using the 0.99 quantile to find reasonable values given a latitude ----
# Q3 + 1.5 IQR gives cutoff values > 1
if(! file.exists('data/avhrr-viirs-ndvi/q99-raster-jan.tif')) {
  for(i in 1:12) {
    print(i)
    month <- tolower(month.abb[i])
    
    r_q99 <-
      map(list.files(
        path = 'data/avhrr-viirs-ndvi/raster-files',
        pattern = if_else(i < 10, paste0('0', i), as.character(i)) %>%
          paste0('.._.*\\.nc'),
        full.names = TRUE) %>%
          `[`(., seq(1, length(.), by = 10)), # take one every 10 days
        \(.fn) .ndvi <- rast(.fn, lyr = 'NDVI')) %>%
      rast() %>%
      quantile(0.99, na.rm = TRUE)
    
    writeRaster(r_q99, paste0('data/avhrr-viirs-ndvi/q99-raster-',
                              month, '.tif'), overwrite = TRUE)
    rm(i, month, r_q99)
  }
}

# plot histograms of 0.99 quantile
png('figures/global-diagnostics/monthly-99th-percentiles.png',
    width = 15, height = 5.5, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12) {
  plot(rast(paste0('data/avhrr-viirs-ndvi/q99-raster-',
                   tolower(month.abb)[i], '.tif')),
       main = month.name[i])
}
dev.off()

# histograms of 0.99 quantile values at high latitudes
# all problematic values are at latitudes > 58
png('figures/global-diagnostics/monthly-99th-percentiles-hist.png',
    width = 15, height = 15, units = 'in', res = 300)
layout(matrix(1:12, ncol = 4, byrow = TRUE))
for(i in 1:12) {
  rast(paste0('data/avhrr-viirs-ndvi/q99-raster-', tolower(month.abb)[i],
              '.tif')) %>%
    as.data.frame(xy = TRUE) %>%
    filter(q0.99 > 0.2, y > 50) %>%
    pull(y) %>%
    hist(main = bquote(bold('Cells with 99'^'th'~'percentile above 0.2 in'~
                              .(month.name[i]))), breaks = 1000,
         xlim = c(50, 85), ylim = c(0, 6e3), xlab = 'Latitude')
  abline(v = 58, col = 'red')
  rm(i)
}
dev.off()

q99 <-
  map(tolower(month.abb), \(.m){
    rast(paste0('data/avhrr-viirs-ndvi/q99-raster-', .m, '.tif')) %>%
      as.data.frame(xy = TRUE) %>%
      mutate(month = .m)
  }) %>%
  bind_rows() %>%
  mutate(month = factor(month, levels = unique(month)))

q99_summary <- q99 %>%
  group_by(month, y) %>%
  summarize(cutoff = quantile(q0.99, 0.1, na.rm = TRUE)) %>%
  mutate(season = case_when(month %in% levels(month)[c(12,1:2)] ~ 'Winter',
                            month %in% levels(month)[c(3:5)] ~ 'Spring',
                            month %in% levels(month)[c(6:8)] ~ 'Summer',
                            month %in% levels(month)[c(9:11)] ~ 'Fall') %>%
           factor(., levels = unique(.)),
         month_int = as.numeric(month),
         month = factor(month.name[month_int], levels = month.name))

# plot 10th percentile of the 90th percentiles against latitude
ggplot() +
  facet_wrap(~ season) +
  geom_vline(xintercept = 60, lty = 'dashed') +
  geom_hline(yintercept = 0.65, lty = 'dashed') +
  annotate(xmin = 60, xmax = Inf, ymin = 0.65, ymax = Inf,
           geom = 'rect', fill = '#FF000030', color = 'red', lwd = 1) +
  geom_line(aes(y, cutoff, color = month, group = month), q99_summary,
            color = 'black', lwd = 1) +
  geom_line(aes(y, cutoff, color = month, group = month), q99_summary) +
  labs(x = 'Latitude',
       y = expression(bold('10'^{th}~'percentile for each latitude of the'~
                             '99'^{th}~'percentiles for each cell'))) +
  scale_color_discreterainbow(name = 'Month') +
  theme(legend.position = 'top')

ggsave('10th-percentile-of-90th-percentiles-across-latitude.png',
       path = 'figures/global-diagnostics', width = 10, height = 10,
       dpi = 300, bg = 'white')

# all problematic values are at latitude > 60
q99_summary %>%
  filter(y > 0, cutoff > 0.7) %>%
  pull(y) %>%
  min()

p_diag_rast <-
  ggplot(q99_summary, aes(month_int, y, fill = cutoff)) +
  geom_raster() +
  geom_hline(yintercept = 60, color = 'red', lty = 'dashed') +
  scale_x_continuous('Month', expand = c(0, 0), breaks = 1:12,
                     labels = month.name) +
  scale_y_continuous('Latitude', expand = c(0, 0), limits = c(-90, 90),
                     breaks = seq(-90, 90, by = 30)) +
  scale_fill_viridis_c(
    expression(bold('10'^{th}~'percentile for each latitude of the'~
                      '99'^{th}~'percentiles for each cell')),
    limits = c(NA, 0.65)) +
  theme(panel.background = element_rect(fill = 'black'),
        legend.position = 'top')

ggsave('10th-percentile-of-90th-percentiles-across-latitude-raster.png',
       plot = p_diag_rast, path = 'figures/global-diagnostics',
       width = 10, height = 10, dpi = 300, bg = 'white')

# find cutoffs for data cleaning ----
# actually, it is best to use time-varying cutoffs of latitude since
# october can still be quite green
p_diag_rast +
  geom_point(aes(color = case_when(month_int == 1 ~ cutoff > 0.2 & y > 60,
                                   month_int == 2 ~ cutoff > 0.2 & y > 60,
                                   month_int == 3 ~ cutoff > 0.2 & y > 60,
                                   month_int == 4 ~ cutoff > 0.2 & y > 60,
                                   month_int == 9 ~ cutoff > 0.2 & y > 70,
                                   month_int == 10 ~ cutoff > 0.2 & y > 70,
                                   month_int == 11 ~ cutoff > 0.2 & y > 60,
                                   month_int == 12 ~ cutoff > 0.2 & y > 60,
                                   .default = FALSE)),
             show.legend = FALSE) +
  scale_color_manual(values = c('transparent', 'red'))

ggsave('10th-percentile-of-90th-percentiles-across-latitude-raster-cutoffs.png',
       path = 'figures/global-diagnostics',
       width = 10, height = 10, dpi = 300, bg = 'white')
