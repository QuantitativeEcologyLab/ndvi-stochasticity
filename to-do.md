# PRIORITY

after predictions are done for neotropic and nearctic:
- cannot predict on lab linux because the NA/SA model alone takes up ~300 GB of RAM
- fit africa, palearctic, then islands, using ~*30 cores* while predicting from other models
- create raster of long-term mean (exclude temporal smooths)
- calculate rasters of observed values, `mu_hat`, and `e^2` for other groups
- fit variance models without `ti(y,x,doy,year)`?
- save daily rasters of estimated mean and DENVar (save for each model and then sum with `na.rm = TRUE`?)
- crete appendix C: global models (code from `010-global-models.R`) and additional figures

# products

rasters of DENVar and corresponding mean estimates:
- daily rasters
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

- update graphical abstract to use fig 1 (Robinson proj) ( https://www.sciencedirect.com/journal/remote-sensing-of-environment/publish/guide-for-authors#toc-27)
- include a concise, descriptive caption for each supplementary file describing its content
- Abbreviate journal names according to the List of Title Word Abbreviations (LTWA)?
- review autor checklist: https://www.sciencedirect.com/journal/remote-sensing-of-environment/publish/guide-for-authors
