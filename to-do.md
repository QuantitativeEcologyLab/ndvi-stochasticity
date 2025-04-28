## tests
- look at green-up rates of areas wit rapid growth rates; e.g.: boreal canada, polar areas
- biomes run diagonally: check that smoooth effects recognize this rather than depending strongly on latitude and produce horizontal contours
- test effects of spatiotemporal aggregation on green-up rates
- test coarser spatial resolution and finer temporal resultion
- run tests for us and canada together (except Hawai'i to since it's so far from the mainland)
- try modeling log(e^2) or using a lognormal distribution if using gaussian model for mean

- check OpenBLAS with openMP on EME linux
- try increasing max dataframe size (ensure `R` is 64-bit): https://stat.ethz.ch/R-manual/R-devel/library/base/html/Memory-limits.html

## modeling
- try running `betals` model on a subset of the data to get an idea of the fitting time
  - test `method = 'NCV'` (faster than REML, but not as fast as `bam()`; see dave miller's work: https://calgary.converged.yt/articles/ncv.html and simon's paper: https://arxiv.org/html/2404.16490v1)
- split by continent to remove issues with monotonicity across greenland with latitude and also increase dataset size
- use `{fmesher}` to construct mesh for MRFs when modeling areas with many islands:
  - https://webhomes.maths.ed.ac.uk/~flindgre/posts/2018-07-22-spatially-varying-mesh-quality/
  - https://wires.onlinelibrary.wiley.com/doi/10.1002/wics.1443
  - use larger triangles near coast to avoid edge effects
  - fit an MRF using vertices of triangles as REs for the spatial smooth
  - use "barriers" to account for separating islands

# products
- static raster for GIS people: exclude `s(year)`, `s(doy)`, and `ti(year,doy)`
- shiny app for predicting from the model given the data
- gif of mean and var over the years and over doy
- integrate shiny app into MoveBank? (use MoveApps?)
- rasters with contours instead of boundaries of ecoregion to check for detection of ecoregions as well as looking for new ones

# figures for paper
- simulated data for mean and var with different smoothness
- plot grid of mean over space, `s(year)`, `s(doy)`, `s(elev)`
- plot grid of DENVar over space, `s(year)`, `s(doy)`, `s(elev)`, `s(mean)`
- histograms of DENVar around flat map of biomes (change colors)
- hex plots of DENVAr with other variables (add temperature, annual precip, extreme weather events, fires, invasive species, max body size, disease occurrence, ...)

supplementary:
- hex plots of pred vs obs for mean and var
- histograms of mean around flat map of biomes (change colors)

