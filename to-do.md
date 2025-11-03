# PRIORITY
- warning for non-matching CRSs is due different CRS for precip layer: remove precip from datasets
```
qread('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-africa-t-4-s-4-ndvi-data.qs', nthreads=62) %>%
  select(! precip_m_day) %>%
  qsave('data/avhrr-viirs-ndvi/group-level-datasets/aggregated-africa-t-4-s-4-ndvi-data.qs',
        nthreads = 62)
```
- finish running mean models
- save rasters of cleaned and aggregated data 
- save rasters of `mu_hat` and `e^2` for each cleaned data raster
- add column of `e^2` to each of the four datasets
- run variance models (use lognormal if normal gives negative values)
- save daily rasters of DENVar
- save rasters of `mu_hat` for days that don't have a data raster

# products

rasters of DENVar and corresponding mean estimates:
- daily rasters
- yearly rasters (exclude `s(doy)`, `ti(year,doy)`, and `ti(y,x,doy)`)
- long-term average rasters (exclude `s(year)`, `s(doy)`, and `ti()` terms)

# data visualization

- gif of mean and var over the years and over doy
- rasters with contour layer instead of boundaries of ecoregions:
  - can we detect the existing biomes?
  - can we identify new biomes based on DENVar?
  - hex plots by group?

# alternative ways of accessing data

- integrate shiny app into MoveBank? (use MoveApps?)
- shiny app extracting values from rasters given a date and coordinates?

# for publication
- add a title in each figure caption
- update graphical abstract to use static rasters of mean NDVI and DENVar ( https://www.sciencedirect.com/journal/remote-sensing-of-environment/publish/guide-for-authors#toc-27)
- update grahpical abstract to use robinson projection in panels B and C
- include a concise, descriptive caption for each supplementary file describing its content
- Abbreviate journal names according to the List of Title Word Abbreviations (LTWA)?
- review autor checklist: https://www.sciencedirect.com/journal/remote-sensing-of-environment/publish/guide-for-authors