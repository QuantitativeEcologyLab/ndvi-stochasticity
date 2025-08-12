rbeta2 <- function(mu, s2) {
  alpha <- ((1 - mu) / s2 - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  
  cat('alpha =', alpha, ', beta =', beta)
  x <- rbeta(1e3, alpha, beta)
  
  hist(x, main = bquote(mu == .(mu)~','~sigma^2 == .(s2)),
       ylab = '', xlab = 'NDVI', xlim = c(0, 1))
}

layout(matrix(1:4, ncol = 2, byrow = TRUE))

# the max value of the variance depends on the mean
rbeta2(mu = 0.5, s2 = 0.01)
rbeta2(mu = 0.5, s2 = 0.25) # scale is 0.25 / (0.5 * 0.5) = 1
rbeta2(mu = 0.9, s2 = 0.02)
rbeta2(mu = 0.9, s2 = 0.089999) # scale is 0.089999 / (0.9 * 0.1) ~= 1

# s2 = 0.01 is approximately the upper limit for reasonable distributions
rbeta2(mu = 0.5, s2 = 0.001)
rbeta2(mu = 0.5, s2 = 0.01)
rbeta2(mu = 0.5, s2 = 0.02)
rbeta2(mu = 0.5, s2 = 0.04)

# s2 > 0.02 easily gives gaussian samples with values < 0 or > 1
hist(rnorm(1e3, mean = 0.5, sd = sqrt(0.001))); abline(v = 0:1, col = 'red')
hist(rnorm(1e3, mean = 0.5, sd = sqrt(0.01))); abline(v = 0:1, col = 'red')
hist(rnorm(1e3, mean = 0.5, sd = sqrt(0.02))); abline(v = 0:1, col = 'red')
hist(rnorm(1e3, mean = 0.5, sd = sqrt(0.04))); abline(v = 0:1, col = 'red')
