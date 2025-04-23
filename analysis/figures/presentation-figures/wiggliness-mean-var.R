library('dplyr') # for data wrangling
library('tidyr') # for data wrangling
library('mgcv') # for modeling
library('ggplot2') # for plotting
source('analysis/figures/000-default-ggplot-theme.R')

d <- tibble(x = runif(1e3),
            mu = sinpi(sqrt(x) * 10 -0.5) + 3 * sqrt(x),
            s2 = cospi(sqrt(x) - 0.5),
            y = rnorm(n = length(x), mean = mu, sd = sqrt(s2)))

layout(1:3); plot(mu ~ x, d); plot(s2 ~ x, d); plot(y ~ x, d); layout(1)

m_1 <- gam(list(y ~ s(x), ~ s(x)), family = gaulss(), data = d, method = 'REML')
l_inv_s <- function(link_pred) (exp(link_pred) + 0.01)

tibble(x = seq(0, 1, by = 1e-3)) %>%
         bind_cols(predict(m_1, newdata = ., type = 'link', se.fit = TRUE) %>%
                     as.data.frame() %>%
                     rename(mu_est = 1, s2_est = 2,
                            mu_se = 3, s2_se = 4)) %>%
  pivot_longer(-x, names_sep = '_', names_to = c('par', 'measure')) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(lab = if_else(par == 'mu', 'Mean', 'Variance'),
         lwr = if_else(par == 'mu', est - 1.96 * se,
                       l_inv_s(est - 1.96 * se)),
         upr = if_else(par == 'mu', est + 1.96 * se,
                       l_inv_s(est + 1.96 * se)),
         est = if_else(par == 'mu', est, l_inv_s(est))) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free', strip.position = 'left', ncol = 1) +
  geom_point(aes(x, y), mutate(d, lab = 'Mean'), alpha = 0.1) +
  geom_ribbon(aes(x, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(x, est), lwd = 1) +
  ylab(NULL) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 12))
  
ggsave('figures/presentation-figures/mean-var-wiggly.png',
       width = 6, height = 6, units = 'in', bg = 'white', dpi = 300)

m_2 <- gam(list(y ~ s(x, k = 4), ~ s(x)), family = gaulss(), data = d, method = 'REML')

tibble(x = seq(0, 1, by = 1e-3)) %>%
  bind_cols(predict(m_2, newdata = ., type = 'link', se.fit = TRUE) %>%
              as.data.frame() %>%
              rename(mu_est = 1, s2_est = 2,
                     mu_se = 3, s2_se = 4)) %>%
  pivot_longer(-x, names_sep = '_', names_to = c('par', 'measure')) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(lab = if_else(par == 'mu', 'Mean', 'Variance'),
         lwr = if_else(par == 'mu', est - 1.96 * se,
                       l_inv_s(est - 1.96 * se)),
         upr = if_else(par == 'mu', est + 1.96 * se,
                       l_inv_s(est + 1.96 * se)),
         est = if_else(par == 'mu', est, l_inv_s(est))) %>%
  ggplot() +
  facet_wrap(~ lab, scales = 'free', strip.position = 'left', ncol = 1) +
  geom_point(aes(x, y), mutate(d, lab = 'Mean'), alpha = 0.1) +
  geom_ribbon(aes(x, ymin = lwr, ymax = upr), alpha = 0.2) +
  geom_line(aes(x, est), lwd = 1) +
  ylab(NULL) +
  theme(strip.placement = 'outside', strip.background = element_blank(),
        strip.text = element_text(size = 12))

ggsave('figures/presentation-figures/mean-var-smooth.png',
       width = 6, height = 6, units = 'in', bg = 'white', dpi = 300)
