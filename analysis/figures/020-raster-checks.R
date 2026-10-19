library('sf')      # for shapefiles
library('terra')   # for rasters
library('dplyr')   # for data wrangling
library('ggplot2') # for fancy plots
library("tidyterra")

save_raster <- function(r, parameter_name, check) {
  r2 <- eval(parse(text = paste("r", check))) %>%
    project(robinson_crs)
  ggplot() +
    tidyterra::geom_spatraster()
}

rm(r, parameter_name, check)

r_mu <- rast("output/long-term-estimates.tif")$mu_hat
r_s2 <- rast("output/long-term-estimates.tif")$s2_hat
r 
# convert to parameters for beta distribution
mu <- (r_mu + 1) / 2 # (-1, 1) to (0, 1)
s2 <- r_s2 / 4

# phi is in (0, 1)
phi <- s2 / (mu * (1 - mu)) # scale parameter
a <- mu * (1 / phi - 1) # shape_1
b <- (1 - mu) * (1 / phi - 1) # shape 2

# locations with impossible shape are all due to negative variance estimates
plot(phi > 0)

layout(matrix(1:4, ncol = 2))
plot(r_mu > -1 & r_mu < 1)
plot(r_s2 < 0)
plot(a < 0)
plot(b < 0)

formatC(sum(values(!is.na(r_mu))), format = "e", digits = 1) # 5.8M cells
any(values(r_mu > 1), na.rm = TRUE) # no cells with NDVI > 1
sum(values(r_mu < -1), na.rm = TRUE) # only 48 cells with NDVI < -1
mean(values(r_mu < -0.1), na.rm = TRUE) # very few with NDVI < -0.1
mean(values((r_s2 < 0) == (b < 0)), na.rm = TRUE) # bad a implies bad s2
mean(values((a < 0) == (b < 0)), na.rm = TRUE) # bad b matches bad b

# find proportion of cells with variance above the possible maximum
# phi = s2 / (mu * (1 - mu))
# s2 = phi * (mu * (1 - mu))
# s2 ranges from 0 to 0.25, so as mu approaches 0.5 phi is maximized at 1
# therefore s2 must be < mu * (1 - mu)
plot(s2 < mu * (1 - mu))
sum(values(s2 > mu * (1 - mu)), na.rm = TRUE) # only 59 cells with impossible variances
values(r_mu)[which(values(s2 > mu * (1 - mu)))] # only 59 cells with impossible variances
