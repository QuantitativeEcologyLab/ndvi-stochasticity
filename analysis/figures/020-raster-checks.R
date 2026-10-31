library('sf')        # for shapefiles
library('terra')     # for rasters
library('dplyr')     # for data wrangling
library('ggplot2')   # for fancy plots
library("tidyterra") # for plotting rasters
source("analysis/figures/000-default-ggplot-theme.R")
source("analysis/figures/000-robinson-objects.R")

options(scipen=10)

save_raster <- function(r, parameter_name, check) {
  r2 <- eval(parse(text = paste("r", check)))
  p <- ggplot() +
    geom_sf(data = bounds, fill = "white") +
    geom_sf(data = shp_rob) +
    geom_spatraster(data = r2, na.rm = TRUE) +
    geom_sf(data = bounds, fill = "transparent", color = "black") +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    khroma::scale_fill_highcontrast(name = paste(parameter_name, check),
                                    na.translate = FALSE) +
    theme(panel.border = element_blank(), legend.position = "top")
  check <- gsub(" ", "", check)
  check <- gsub("<", " lt ", check)
  check <- gsub(">", " gt ", check)
  ggsave(paste0(parameter_name, check, ".png"), p,
         path = "figures/global-diagnostics/appendix-c",
         width = 6, height = 4, dpi = 300, bg = "white")
}

r_mu <- rast("output/long-term-estimates.tif")$mu_hat
r_s2 <- rast("output/long-term-estimates.tif")$s2_hat

# convert to parameters for beta distribution
mu <- (r_mu + 1) / 2 # (-1, 1) to (0, 1)
s2 <- r_s2 / 4

# phi is in (0, 1)
phi <- s2 / (mu * (1 - mu)) # scale parameter
a <- mu * (1 / phi - 1) # shape_1
b <- (1 - mu) * (1 / phi - 1) # shape 2

# summaries for unprojected rasters ----
cli::cli_alert_info(paste0(
  "The unprojected rasters of long-term mean NDVI and DENVar had a total of ",
  round(unique(c(sum(values(!is.na(r_mu))),
                 sum(values(!is.na(r_mu))))) / 1e6, 1), " million non-NA ",
  "cells each."))

cli::cli_alert_info(paste0(
  "The raster of mean NDVI had ",
  sum(values(r_mu > 1), na.rm = TRUE), " cells with values > 1 and ",
  sum(values(r_mu < -0.1), na.rm = TRUE), " cells (",
  round(mean(values(r_mu < -0.1), na.rm = TRUE) * 100, 3),
  "%) with NDVI < -0.1. Only ",
  sum(values(r_mu < -1), na.rm = TRUE), " cells (",
  round(mean(values(r_mu < -1), na.rm = TRUE) * 100, 4),
  "%) had NDVI < -1."))

stopifnot(all(values((r_mu < -1) == (r_s2 > 1)), na.rm = TRUE))
cli::cli_alert_info(paste0(
  "The DENVar raster had ",
  sum(values(r_s2 > 1), na.rm = TRUE), " cells with values > 1 and ",
  sum(values(r_s2 < 0), na.rm = TRUE), " cells (",
  round(mean(values(r_s2 < 0), na.rm = TRUE) * 100, 3),
  "%) with values < 0. All cells with DENVar values > 1 had mean NDVI < -1."))

cli::cli_alert_info(paste0(
  "After scaling the rasters to be within the [0, 1] range for NDVI and ",
  "(0, 0.25) range for DENVar, and calculating the corresponding ",
  "parameters for a beta distribution (\U03B1 > 0, \U03B2 > 0, \U03C6 ",
  "\U2208 (0, 1)): ",
  formatC(sum(values(a < 0), na.rm = TRUE), big.mark = ","), " cells (",
  round(mean(values(a < 0), na.rm = TRUE) * 100, 3),
  "%) had negative estimated shape paramter \U03B1. ",
  formatC(sum(values(b < 0), na.rm = TRUE), big.mark = ","), " cells (",
  round(mean(values(b < 0), na.rm = TRUE) * 100, 3),
  "%) had negative estimated shape paramter \U03B2. ",
  formatC(sum(values(phi < 0), na.rm = TRUE), big.mark = ","), " cells (",
  round(mean(values(phi < 0), na.rm = TRUE) * 100, 3),
  "%) had negative estimated scale paramter \U03C6, and ",
  formatC(sum(values(phi < 0), na.rm = TRUE), big.mark = ","), " cells (",
  round(mean(values(phi < 0), na.rm = TRUE) * 100, 3),
  "%) had an estimated scale paramter \U03C6 > 1. ",
  round(mean(values((a < 0) == (b < 0)), na.rm = TRUE) * 100, 4),
  "% of cells with negative \U03B1 also had a negative \U03B2."
))

# figures of rasters projected to Robinson CRS ----
r_mu <- r_mu %>% project(robinson_crs)
r_s2 <- r_s2 %>% project(robinson_crs)
mu <- mu %>% project(robinson_crs)
s2 <- s2 %>% project(robinson_crs)
phi <- phi %>% project(robinson_crs)
a <- a %>% project(robinson_crs)
b <- b %>% project(robinson_crs)

# locations with impossible shape are all due to negative variance estimates
save_raster(r_mu, "Mean NDVI", "< -1")
save_raster(r_mu, "Mean NDVI", "> 1")
save_raster(r_s2, "DENVar", "< 0")
save_raster(mu, "Scaled mean NDVI", "< 0")
save_raster(mu, "Scaled mean NDVI", "> 1")
save_raster(s2, "Scaled DENVar", "< 0")
save_raster(s2, "Scaled DENVar", "> 0.25")
save_raster(a, "\U03B1", "< 0")
save_raster(b, "\U03B2", "< 0")
save_raster(phi, "\U03C6", "< 0")
save_raster(phi, "\U03C6", "> 1")
