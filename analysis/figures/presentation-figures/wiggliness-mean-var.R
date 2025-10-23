#' caption (*NOTE: MISSING METHODS*)
"Estimates of variance in NDVI depend on the estimated trends in mean NDVI.} Black lines indicate the estimated mean and variance in simulated NDVI (points). Grey ribbons correspond to 95% credible intervals, while the colored segments highlight significant increases (red) and decreases (blue).\\\\ Smoother estimates of the mean (top left) result in wigglier and larger variance estimates (bottom left), since short-term trends in the data are considered to be stochastic events. In contrast, wigglier estimates of the mean (top right) result in smoother and smaller variance estimates (bottom right), since short-term trends in the data are attributed to changes in the mean. The data in the top row is the same in both columns. It was generated as the sum of a sinusoidal function and a mean-reverting, correlated, stochastic process with mean zero."

library('dplyr')   # for data wrangling
library('tidyr')   # for data wrangling
library('mgcv')    # for modeling
library('gratia')  # for working with GAMs
library('ggplot2') # for plotting
source('analysis/figures/000-default-ggplot-theme.R')

# OU process given samples times, reversion time, standard deviation 
ou <- function(t, theta, sigma){
  n <- length(t) # number of samples
  dw  <- rnorm(n = n, mean = 0, sd = sqrt(sigma)) # noise between steps
  
  x <- numeric(n)
  x[1] <- 0 # set starting point
  for (i in 2:n) {
    x[i] <- (1 - theta) * x[i-1] + dw[i - 1]
  }
  return(x)
}

set.seed(3)
d <- tibble(t = seq(0, 1, length.out = 500),
            mu = (sinpi(sqrt(t + 0.5) * 10 - 0.5) * 3 + 1 * sqrt(t)) / 50 + 0.3,
            s2 = cospi(sqrt(t) - 0.5) * 0.001^2,
            noise = ou(t = t, theta = 0.2, sigma = sqrt(s2)),
            y = mu + noise)

if(FALSE) {
  ggplot(d) +
    geom_line(aes(t, mu), color = 'red') +
    geom_line(aes(t, y), alpha = 0.5)
}

m_wiggly <- gam(list(y ~ s(t, k = 30, bs = 'ad'), ~ s(t)),
                family = gaulss(), data = d, method = 'REML')
m_smooth <- gam(list(y ~ s(t, k = 10), ~ s(t)),
                family = gaulss(), data = d, method = 'REML')

get_preds <- function(model, model_name) {
  d0 <- tibble(t = seq(0, 1, by = 1e-3))
  
  inv_link_l <- inv_link(model, parameter = 'location')
  inv_link_s <- inv_link(model, parameter = 'scale')
  
  bind_cols(
    d0,
    predict(model, newdata = d0, type = 'link', se.fit = TRUE) %>%
      as.data.frame() %>%
      transmute(Mean_est = inv_link_l(fit.1),
                Mean_lwr = inv_link_l(fit.1 - 1.96 * se.fit.1),
                Mean_upr = inv_link_l(fit.1 + 1.96 * se.fit.1),
                Variance_est = 1 / inv_link_s(fit.2)^2, # (SD)^2
                Variance_lwr = 1 / inv_link_s(fit.2 - 1.96 * se.fit.2)^2,
                Variance_upr = 1 / inv_link_s(fit.2 + 1.96 * se.fit.2)^2,
                model = model_name)) %>%
    pivot_longer(-c(t, model), names_sep = '_',
                 names_to = c('parameter', 'measure')) %>%
    mutate(parameter = if_else(parameter == 'Mean', 'Mean NDVI',
                               'Variance in NDVI')) %>%
    pivot_wider(names_from = 'measure', values_from = 'value') %>%
    arrange(parameter, t) %>% # so mean and var match with s(t) and s.1(t)
    bind_cols(.,
              derivatives(model, select = c('s(t)', 's.1(t)'), data = d0) %>%
                transmute(smooth = .smooth,
                          deriv = .derivative,
                          deriv_lwr = .lower_ci,
                          deriv_upr = .upper_ci,
                          trend = case_when(deriv_lwr > 0 ~ 'increase',
                                            deriv_upr < 0 ~ 'decrease')))
}

preds <-
  bind_rows(get_preds(m_wiggly, 'Wiggly'), get_preds(m_smooth, 'Smooth')) %>%
  mutate(model = paste(model, 'mean'),
         # to keep highlights distinct between changes in direction
         g = c(0, cumsum(trend[1:(n()-1)] != trend[2:n()])))

resids <- bind_rows(
  mutate(d,
         model = 'Smooth mean',
         e = residuals(m_smooth, type = 'response')),
  mutate(d,
         model = 'Wiggly mean',
         e = residuals(m_wiggly, type = 'response'))) %>%
  mutate(parameter = 'Variance in NDVI')

fig <-
  ggplot(preds) +
  facet_grid(parameter ~ model, scales = 'free', switch = 'y', ) +
  # data and squared residuals
  geom_point(aes(t, y), mutate(d, parameter = 'Mean NDVI'), alpha = 0.2) +
  geom_point(aes(t, e^2), resids, alpha = 0.2) +
  # estimated trends
  geom_ribbon(aes(t, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(t, est, color = trend, group = g), linewidth = 2.5,
            alpha = 1, show.legend = FALSE, na.rm = TRUE) +
  geom_line(aes(t, est), linewidth = 0.75) +
  ylab(NULL) +
  scale_x_continuous('Time', labels = NULL) +
  scale_color_manual(values = c('#CA0020', '#2166AC'), na.value = NA,
                     breaks = c('increase', 'decrease')) +
  theme(strip.placement = 'outside', strip.background.y = element_blank(),
        strip.text = element_text(size = 12),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12))

TEST_CDV <- FALSE
if(TEST_CDV) {
  colorblindr::cvd_grid(fig) # test for color-vision deficiency
} else {
  fig
}

ggsave('figures/presentation-figures/mean-var-smooth-wiggly.png', fig,
       width = 6, height = 6, units = 'in', bg = 'white', dpi = 300)
