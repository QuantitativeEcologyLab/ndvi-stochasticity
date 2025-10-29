# PRIORITY
- warning for non-matching CRSs is due different CRS for precip layer: remove precip from datasets
```
qread('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-africa-t-4-s-4-ndvi-data.qs', nthreads=62) %>%
  select(! precip_m_day) %>%
  qsave('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-africa-t-4-s-4-ndvi-data.qs',
        nthreads = 62)
```
- save the aggregated and cleaned rasters
- run: (1) americas, (2) afrotropics, (3) islands, (4) paleoarctic

# tests
- model e^2 using a lognormal distribution

# products
rasters:
- daily rasters
- yearly rasters (exclude `s(doy)`, `ti(year,doy)`, and `ti(y,x,doy)`)
- 2 long-term average rasters: exclude `s(year)`, `s(doy)`, and `ti()` terms 

viz:
- gif of mean and var over the years and over doy
- rasters with contour layer instead of boundaries of ecoregions:
  - can we detect the existing biomes?
  - can we identify new biomes based on DENVar?
  - hex plots by group?

alternative ways of accessing data:
- integrate shiny app into MoveBank? (use MoveApps?)
- shiny app extracting values from rasters?
