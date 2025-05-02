library('mgcv')    # for modeling
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
source('analysis/figures/000-default-ggplot-theme.R')

d <- tidyr::expand_grid(cell_id = 1:4,
                        time = 1:1e3) %>%
  mutate(mu = sin(time / 30) * sqrt(cell_id),
         fact_id = factor(cell_id),
         y = rnorm(n(), mean = mu, sd = 0.2))
d

extra <- list(nbs = list(`1` = c(0),
                         `2` = c(3, 4),
                         `3` = c(2, 4),
                         `4` = c(2, 3)))

ggplot(d) +
  facet_wrap(~ fact_id) +
  geom_point(aes(time, y), alpha = 0.1) +
  geom_line(aes(time, mu), color = 'darkorange', lwd = 1)

m0 <- bam(y ~ fact_id + s(time, by = fact_id, k = 30),
          data = d, method = 'fREML', discrete = TRUE)
gratia::draw(m0, rug = FALSE, n = 1e3)

m1 <-
  bam(y ~ s(time, k = 30) +
        s(fact_id, k = 4, bs = 'mrf', xt = extra) +
        te(time, fact_id, k = c(10, 4), bs = c('tp', 'mrf'), xt = extra),
      data = d, method = 'fREML', discrete = TRUE)

m2 <-
  bam(y ~ te(time, fact_id, k = c(30, 4), bs = c('tp', 'mrf'), xt = extra),
      data = d, method = 'fREML', discrete = TRUE)

d %>%
  mutate(by_smooth = predict(m0),
         `s(mrf) and te(mrf)` = predict(m1),
         `te(mrf) only` = predict(m2)) %>%
  tidyr::pivot_longer(cols = by_smooth:`te(mrf) only`) %>%
  ggplot() +
  facet_grid(name ~ fact_id) +
  geom_point(aes(time, y), alpha = 0.1) +
  geom_line(aes(time, value, color = name), lwd = 1) +
  khroma::scale_color_highcontrast() +
  theme(legend.position = 'none')
