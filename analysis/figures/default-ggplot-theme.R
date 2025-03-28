library('ggplot2') # for fancy figures
library('khroma')  # for colorblind-friendly palettes

.PLOT_PALETTES <- FALSE

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(face = 'bold')))

# custom NDVI color palette
ndvi_pal <- plot_scheme(color('bukavu')(30)) # wrong order
if(.PLOT_PALETTES) plot_scheme(ndvi_pal)
ndvi_pal <- ndvi_pal[c(1:13, length(ndvi_pal):(length(ndvi_pal)/2 + 1))]
if(.PLOT_PALETTES) plot_scheme(ndvi_pal)
create_ndvi_pal <- colorRampPalette(ndvi_pal)
ndvi_pal <- create_ndvi_pal(100)
if(.PLOT_PALETTES) plot_scheme_colorblind(ndvi_pal)
if(.PLOT_PALETTES) {
  expand.grid(x = 1:10,
              y = 1:10) %>%
    mutate(z = runif(100, min = -1, max = 1)) %>%
    ggplot() +
    geom_raster(aes(x, y, fill = z)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradientn(colours = ndvi_pal, limits = c(-1, 1))
  
  (expand.grid(x = 1:10,
               y = 1:10) %>%
     mutate(z = runif(100, min = -1, max = 1)) %>%
     ggplot() +
     geom_raster(aes(x, y, fill = z)) +
     scale_x_continuous(expand = c(0, 0)) +
     scale_y_continuous(expand = c(0, 0)) +
     scale_fill_gradientn(colours = ndvi_pal, limits = c(-1, 1))) %>%
    colorblindr::cvd_grid()
}

