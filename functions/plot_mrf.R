library('mgcv')
library('ggplot2')
library('dplyr')

plot_mrf <- function(.model, .term, .newdata, .type = 'response',
                     .fun = ndvi_to_11, .limits = c(-1, 1)) {
  # find unique cells
  .newdata <- .newdata %>%
    group_by(x, y) %>%
    slice(1) %>%
    ungroup()
  
  # find predictions
  .newdata <- .newdata %>%
    mutate(mu_hat = predict(.model, newdata = .newdata, terms = .term,
                            type  = .type, se.fit = FALSE) %>%
             .fun())
  
  # apply transformation, if necessary
  ggplot(.newdata, aes(x, y, fill = mu_hat)) +
    coord_equal() +
    geom_raster() +
    labs(x = NULL, y = NULL) +
    scale_fill_gradientn('NDVI', colours = ndvi_pal, limits = .limits)
}
