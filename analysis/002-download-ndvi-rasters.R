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
              'data/avhrr-viirs-ndvi/raster-files/')
if(! dir.exists(DIR)) stop('DIR does not exist!')
plan() # check plan

output <-
  future_map2_chr(fn_tib$years, fn_tib$fn, function(.year, .filename) {
    if(! file.exists(paste0(DIR, .filename))) {
      download.file(
        url = paste0(url_main, .year, '/', .filename),
        #' cannot save directly to `//home/shared` because it needs sudo access
        destfile = paste0(DIR, .filename), method = 'libcurl', quiet = TRUE)
      return(paste('Downloaded', .filename))
    } else {
      warning(paste('Did not download', .filename))
      return(paste('Did not download', .filename))
    }
  }, .progress = TRUE)
table(grepl('Downloaded', output))

#' There's a few days that have no raster associated with them, but it's
#' generally rare. We have 8 2-day gaps, 3 3-day gaps, 2 4-day gaps, and
#' one 16-day gap. The gaps range from 1982 to 2024, and the longest one
#' is in summer of 2022.

# check files
fn_dates <- fn_tib %>%
  pull(fn) %>%
  substr(.,
         nchar(.) - nchar('YYYYMMDD_c20240123205954.nc') + 1,
         nchar(.) - nchar('_c20240123205954.nc')) %>%
  as.Date(format = '%Y%m%d')

main_dates <-
  list.files('data/avhrr-viirs-ndvi/raster-files', pattern = '*.nc',
             recursive = FALSE) %>%
  substr(.,
         nchar(.) - nchar('YYYYMMDD_cYYYYmmddHHMMSS.nc') + 1,
         nchar(.) - nchar('_cYYYYmmddHHMMSS.nc')) %>%
  as.Date(format = '%Y%m%d')

# ensure date ranges are the same
range(fn_dates)
all(range(fn_dates) == range(main_dates))

# ensure all rasters downloaded
sum(! fn_dates %in% main_dates)

# check that all files were downloaded
gaps <-
  tibble(date = sort(main_dates),
         before = lag(date),
         gap = date - before,
         year = lubridate::year(date)) %>%
  relocate(before, .before = 1)
table(gaps$gap)
gaps %>%
  filter(gap > 1) %>%
  View()

length(fn_dates) == length(main_dates)
length(fn_dates) - length(main_dates)
as.numeric(diff(range(main_dates), units = 'days')) / 365

#' *need to move files to "//home/shared" folder on the lab Linux*
