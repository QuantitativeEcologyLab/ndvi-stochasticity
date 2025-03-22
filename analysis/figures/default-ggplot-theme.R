library('ggplot2') # for fancy figures

theme_set(theme_bw() +
            theme(panel.grid = element_blank(),
                  text = element_text(face = 'bold')))


# custom NDVI color palette
create_ndvi_pal <- colorRampPalette(c('darkblue', 'dodgerblue', '#744700',
                                      '#d9bb94', 'darkgreen'))
ndvi_pal <- create_ndvi_pal(100)

if(FALSE) khroma::plot_scheme_colorblind(ndvi_pal)
