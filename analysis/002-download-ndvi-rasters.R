library('rvest') # for harvesting web pages; see https://rvest.tidyverse.org/
library('ncdf4') # for nc rasters
library('dplyr') # for data wrangling
library('purrr') # for functional programming
library('furrr') # for parallel functional programming
ncores <- min(50, future::availableCores(logical = FALSE) / 2)
plan(multisession, workers = ncores)

# main url and years to find all folders in the parent directory
url_main <- 'https://www.ncei.noaa.gov/data/land-normalized-difference-vegetation-index/access/'
years <- 1981:2025 #' years to download; *downloading up to 2025-03-16*

#' *\/ delete*
# rasters that Rekha had downloaded that had "preliminary" in the file name
# need to re-download
prelim <- list.files('//./home/shared/NOAA_Files/preliminary-rasters/',
                     include.dirs = FALSE) %>%
  `[`(., which(grepl('preliminary', .))) %>%
  gsub('VIIRS-Land_v001-preliminary_NPP13C1_S-NPP_', '', .) %>%
  gsub('AVHRR-Land_v005-preliminary_AVH13C1_NOAA-19_', '', .) %>%
  substr(., 1, nchar('yyyymmdd')) %>%
  as.Date(format = '%Y%m%d')

length(prelim)
range(prelim)
years <- 2014:2025 #' re-downloading files that were preliminary
#' */\ delete*

# find file names for each raster
fn_tib <- tibble(
  years = years,
  urls = future_map(years, function(y) {
    url_y <- paste0(url_main, y, '/')
    daysfiles <- XML::htmlParse(httr::GET(url_y)) %>%
      XML::xpathSApply(., path = '//a', XML::xmlGetAttr, 'href')
    
    # identify and subset all of the links with '.nc' in the name
    return(tibble(fn = daysfiles[grepl('\\.nc', daysfiles)]))
  })) %>%
  tidyr::unnest(urls)
fn_tib
tail(fn_tib, 10)
sum(grepl('preliminary', fn_tib$fn)) # check if any rasters are preliminary

# download the files for each year
options('timeout') # in seconds
options(timeout = 60 * 10) # increase timeout to 10 minutes
options('timeout') # in seconds

DIR <- paste0('//home/mezzinis/H/GitHub/ndvi-stochasticity/',
              'data/avhrr-viirs-ndvi/raster-files/new/')
if(! dir.exists(DIR)) stop('DIR does not exist!')
plan() # check plan

output <-
  future_map2_chr(fn_tib$years, fn_tib$fn, function(.year, .filename) {
    if(! file.exists(paste0(DIR, .filename))) {
      download.file(
        url = paste0(url_main, .year, '/', .filename),
        # cannot save directly to '//home/shared' because it needs sudo access
        destfile = paste0(DIR, .filename), method = 'libcurl', quiet = TRUE)
      return(paste('Downloaded', .filename))
    } else {
      warning(paste('Did not download', .filename))
      return(paste('Did not download', .filename))
    }
  }, .progress = TRUE)
table(grepl('Downloaded', output))

# move files to the main folder (out of "/new")
new <- list.files(DIR)
length(new)
file.copy(from = paste0(DIR, new), to = paste0(gsub('new/', '', DIR), new))

# check files
new_dates <- list.files('data/avhrr-viirs-ndvi/raster-files/new') %>%
  substr(.,
         nchar(.) - nchar('YYYYMMDD_c20240123205954.nc') + 1,
         nchar(.) - nchar('_c20240123205954.nc')) %>%
  as.Date(format = '%Y%m%d')

main_dates <-
  list.files('data/avhrr-viirs-ndvi/raster-files', pattern = '*.nc') %>%
  substr(.,
         nchar(.) - nchar('YYYYMMDD_c20240123205954.nc') + 1,
         nchar(.) - nchar('_c20240123205954.nc')) %>%
  as.Date(format = '%Y%m%d')

range(new_dates)
range(main_dates)

#' move files to "//home/shared" folder  on the lab Linux *need to do*
new <- list.files(DIR)
length(new)
file.copy(from = paste0(DIR, new), to = paste0(gsub('new/', '', DIR), new))
