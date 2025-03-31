rbeta2 <- function(mu, s2) {
  alpha <- ((1 - mu) / s2 - 1/mu) * mu^2
  beta <- alpha * (1/mu - 1)
  
  cat('alpha =', alpha, ', beta =', beta)
  x <- rbeta(1e3, alpha, beta)
  
  hist(x, main = bquote(mu == .(mu)~','~sigma^2 == .(s2)),
       ylab = '', xlab = 'NDVI', xlim = c(0, 1))
}

layout(matrix(1:4, ncol = 2))
rbeta2(mu = 0.5, s2 = 0.01)
rbeta2(mu = 0.5, s2 = 0.25) # note: scale is 0.25 / (0.5 * 0.5) = 1
rbeta2(mu = 0.9, s2 = 0.02)
rbeta2(mu = 0.9, s2 = 0.089999) # note: scale is 0.089999 / (0.9 * 0.1) = 0.9999889 ~= 1

rbeta2(mu = 0.5, s2 = 0.002)
rbeta2(mu = 0.5, s2 = 0.02) # note: scale is 0.25 / (0.5 * 0.5) = 1
rbeta2(mu = 0.9, s2 = 0.02)
rbeta2(mu = 0.1, s2 = 0.02) # note: scale is 0.089999 / (0.9 * 0.1) = 0.9999889 ~= 1
