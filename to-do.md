# PRIORITY
- run afrotropics next, then islands
- save and import rds files using `qs::qsave()` and `qs::qread()`

# tests
- model e^2 using a lognormal distribution

# products
- static raster for GIS people: exclude `s(year)`, `s(doy)`, and `ti(year,doy)`
- shiny app for predicting from the model given the data
- gif of mean and var over the years and over doy
- integrate shiny app into MoveBank? (use MoveApps?)
- rasters with contours instead of boundaries of ecoregion:
  - check for detection of biomes
  - look for new biomes based on DENVar
  - hex plots by group?
