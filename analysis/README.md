This folder contains all `R` scripts used in this project. The scripts are numbered in the order they should be ran in. The scripts in `analysis/figures` follow the same numbering as the scripts in `anlysis`.

* `00X`: setup
  * `001-download-ndvi-rasters.R` downloads all the NDVI data
  * `002-sardinia-test.R` runs some tests for model setup using the Mediterranean island of Sardigna (Sardinia)
  * `003-ecoregion-polygons.R` creates a shapefile of ecoregion polygons used throughout the rest of the analyses
  * `004-arctic-test.R` runs some tests using an island in Inuit Nunangat. This script illustrates the short growth seasons of arctic regions and some extremely high NDVI values in winter 
  * `005-check-high-lat-data.R` runs some tests for removing excessively high NDVI values at high latitudes in winter
  * `006-taiga-test.R` runs some tests using a large taiga region that was historically and still is inhabited by many different Indigenous Peoples. This script illustrates that the `mrf` smooth is not feasible for large datasets
  * `007-merging-rasters.R` merges the NDVI rasters into single-dataframe datasets for each study region
  * `008-cell-neighbors.R` creates a list of each cell's neighboring cells, which is used for markov random fields smooths
  * `009-test-hgams.R` fits some preliminary test models on the global dataset
* `02X`: modeling mean NDVI
  * `011-mean-ndvi-hgam.R`
* `03X`: modeling variance in NDVI for a given mean
  * `012-var-ndvi-hgam.R`
* other folders
  * `figures` contains `R` scripts that are specifically for creating figures
  * `misc` contains some miscellaneous `R` scripts for minor tests
