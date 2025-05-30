# this function decodes the signed short (16-bit with +/-) used for quality
# assessment in the QA layer of each raster
#' for additional info and examples, see:
#'  - `https://github.com/rspatial/luna/blob/f7bc1cb47dde54a75992b1c74ca6e7f7ddd2a9fc/R/modis_qc.R#L4`
#'  - `https://github.com/rspatial/luna/blob/f7bc1cb47dde54a75992b1c74ca6e7f7ddd2a9fc/man/modis_mask.Rd`
library('dplyr')

decode_qa <- function(qa_values, return_bit = FALSE, warn = FALSE) {
  all_qa_values <- qa_values
  qa_values <- unique(all_qa_values)
  
  if(return_bit) {
    out <-
      furrr::future_map(qa_values, function(.qa) {
        #' to calculate a value from bits, see `functions/bits_to_int.R`
        return(tibble(qa_values = .qa,
                      true_bit = as.numeric(intToBits(.qa)[16:1]) %>%
                        paste(collapse = ''),
                      binary = as.numeric(true_bit),
                      R_bit = as.numeric(intToBits(.qa)[1:16]) %>%
                        paste(collapse = '')))
      }) %>%
      bind_rows()
  } else {
    out <-
      furrr::future_map(qa_values, function(.qa) {
        if(! is.finite(.qa)) return(data.frame(qa = .qa, bit = NA))
        
        # extract a logical version of the short (indices works L -> R)
        k <- as.logical(intToBits(.qa)[1:16])
        
        unused <- c(7, 11:14) + 1 # bits start at 0, but indices start at 1
        
        if(warn){
          if(any(k[unused])) {
            warning(
              paste0(
                'Ignored ', sum(k[unused]), ' unused flags (',
                # list unused flags that were equal to 1 (count from 0)
                paste(unused[which(k[unused])] - 1, collapse = ', '),
                ') that had been set to "1" in short ',
                paste(as.numeric(k), collapse = ''),
                ' = ', .qa, '.')
            )
          }
        }
        
        #' extract meaning from `k`: bits start at 0, but indices start at 1
        #' note: there are `4*2*5*2*2*2*2*2 = 1280` possible numbers
        data.frame(
          qa = .qa,
          bit = paste(as.numeric(k), collapse = ''),
          #' 0-1 Cloud State
          #'    00 Confident Clear
          #'    01 Probably Clear
          #'    10 Probably Cloudy
          #'    11 Confident Cloudy
          cloud_state = case_when(! k[1] & ! k[2] ~ 'Confident Clear',
                                  k[1] ~ 'Probably Clear',
                                  k[2] ~ 'Probably Cloudy',
                                  k[1] & k[2] ~ 'Confident Cloudy'),
          #' 2 Cloud shadow
          #'    1 Yes
          #'    0 No
          cloud_shadow = if_else(k[3], true = 'Yes', false = 'No'),
          #' 3-5 Land/Water
          #'    000 Land & Desert
          #'    001 Land no desert
          #'    010 Inland Water
          #'    011 Sea Water
          #'    100 ---
          #'    101 Coastal
          #'    110 ---
          #'    111 ---
          land_type = case_when(
            ! k[4] & ! k[5] & ! k[6] ~ 'Land & Desert',
            ! k[4] & ! k[5] &   k[6] ~ 'Land no desert',
            ! k[4] &   k[5] & ! k[6] ~ 'Inland Water',
            ! k[4] &   k[5] &   k[6] ~ 'Sea Water',
            k[4] & ! k[5] & ! k[6] ~ NA_character_,
            k[4] & ! k[5] &   k[6] ~ 'Coastal',
            k[4] &   k[5] & ! k[6] ~ NA_character_,
            k[4] &   k[5] &   k[6] ~ NA_character_),
          #' 6 Overall Aerosol Quality
          #'    1 OK
          #'    0 Poor
          aerosol_quality = if_else(k[7], true = 'OK', false = 'Poor'),
          #' 7 Unused
          #' 8 Thin cirrus reflective
          #'    1 Yes
          #'    0 No
          thin_cirrus_reflective = if_else(k[9], true = 'Yes', false = 'No'),
          #' 9 Thin cirrus emissive
          #'    1 Yes
          #'    0 No
          thin_cirrus_emissive = if_else(k[10], true = 'Yes', false = 'No'),
          #' 10 Cloud flag
          #'    1 Cloud
          #'    0 No cloud
          cloud_flag = if_else(k[11], true = 'Cloud', false = 'No cloud'),
          #' 11-14 Unused
          #' 15 Snow/Ice Flag
          #'    1 Snow/Ice
          #'    0 No snow/Ice
          snow_ice = if_else(k[16], true = 'Snow/Ice', false = 'No snow/Ice'))
      }) %>%
      bind_rows()
    
    out <- tibble(qa = all_qa_values) %>%
      left_join(out, by = 'qa')
  } # close else
  
  return(out)
}

if(FALSE) {
  library('terra')
  library('sf')
  source('functions/bit_to_int.R')
  
  # check that values can be back-transformed
  strsplit(decode_qa(10, return_bit = TRUE)$R_bit[1], split = '*')[[1]] %>%
    as.numeric() %>%
    bit_to_int()
  
  decode_qa(c(0:6), return_bit = TRUE) # check shorts
  decode_qa(0) # example of a possible flag
  decode_qa(bit_to_int(c(rep(0, 7), 1, rep(0, 8)))) # with an unused flag
  c(rep(0, 7), 1, rep(0, 2), 1, 1, 1, 1, 0, 0) %>% # multiple unused flags
    bit_to_int() %>%
    decode_qa()
  decode_qa(c(32, 33)) # example of multiple possible flags
  
  # use an exaple from a raster
  sahara <- st_read('data/ecoregions/ecoregions-polygons.shp') %>%
    filter(ECO_NAME == 'Sahara Desert') %>%
    st_geometry() %>%
    st_as_sf()
  
  r <- rast('data/avhrr-viirs-ndvi/raster-files/VIIRS-Land_v001_JP113C1_NOAA-20_20250507_c20250513122857.nc') %>%
    crop(sahara, mask = TRUE)
  plot(r)
  
  decoded <- decode_qa(as.numeric(values(r$QA)))
  
  # make rasters for each QA parameter
  for(cn in colnames(decoded)[3:ncol(decoded)]) {
    r[[cn]] <- r$QA
    values(r[[cn]]) <- decoded[[cn]]
    r[[cn]] <- mask(r[[cn]], r$QA)
  }
  plot(r)
  
  # the 7th flag is occasionally either 1 or 0 despite not being used...
  map(strsplit(filter(decoded, ! is.na(qa))$bit, '*'), \(l) l[7 + 1]) %>%
    unique()
  
  # and it's not due to incorrect reversal of the bit
  map(strsplit(filter(decoded, ! is.na(qa))$bit, '*'), \(l) rev(l)[7 + 1]) %>%
    unique()
}
