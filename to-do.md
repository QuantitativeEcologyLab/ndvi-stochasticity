# PRIORITY

- test aggregation with sardinia example (plot preds vs preds for original model and aggregated):
  - DS
  - MRF
  - MRF aggr
  - bottom row:
    - comparison of preds vs preds (not aggr vs aggr)
    - raster of difference of spatial terms
    - var estimates of two models

- run afrotropics next, then islands

- need to re-run all scripts:
  - create objects of global cell ids and neighbors 
  - split data by continents (not splitting neighbors by continents)
  - re-organize and clean up scripts

# tests

- north america test:
  - create canada data and test for single-day mrf
  - check that mrfs are giving good estimates for biomes and coastlines
  - more complex models
  - look at green-up rates of areas with rapid growth rates; e.g.: boreal canada, polar areas
  - biomes run diagonally: check that smoooth effects recognize this rather than depending strongly on latitude and produce horizontal contours
  - test effects of spatiotemporal aggregation on green-up rates
  - try modeling e^2 using a lognormal distribution
  - test coarser spatial resolution and finer temporal resultion
- try increasing max dataframe size (ensure `R` is 64-bit): https://stat.ethz.ch/R-manual/R-devel/library/base/html/Memory-limits.html

# modeling

- vector size on EME linux is not limited to `2^32 - 1`
- split by continent to remove issues with monotonicity across greenland with latitude and also increase dataset size
- use adaptive splines for year to allow for rapid change in ~2010? (no difference in sardinia test)

# products

- static raster for GIS people: exclude `s(year)`, `s(doy)`, and `ti(year,doy)`
- shiny app for predicting from the model given the data
- gif of mean and var over the years and over doy
- integrate shiny app into MoveBank? (use MoveApps?)
- rasters with contours instead of boundaries of ecoregion:
  - check for detection of biomes
  - look for new biomes based on DENVar
  - hex plots by group?
