# PRIORITY

- split data by continents
- run afrotropics next, then islands

# tests

- model e^2 using a lognormal distribution
- try increasing max dataframe size (ensure `R` is 64-bit): https://stat.ethz.ch/R-manual/R-devel/library/base/html/Memory-limits.html

# products

- static raster for GIS people: exclude `s(year)`, `s(doy)`, and `ti(year,doy)`
- shiny app for predicting from the model given the data
- gif of mean and var over the years and over doy
- integrate shiny app into MoveBank? (use MoveApps?)
- rasters with contours instead of boundaries of ecoregion:
  - check for detection of biomes
  - look for new biomes based on DENVar
  - hex plots by group?
