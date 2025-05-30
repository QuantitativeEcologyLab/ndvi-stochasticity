# a function to convert bits to integers/numeric
#' *NOTE: this function uses R-type (i.e., reversed) bits!*
bit_to_int <- function(bits) {
  if(class(bits) == 'raw') {
    bits <- as.numeric(bits)
  } else if(! class(bits) %in% c('numeric', 'integer')) {
    stop('`bits` must be of class "raw", "numeric" or "integer".\n')
  }
  as.numeric(t(2^(0:14)) %*% bits[1:15] - 2^15 * bits[16])
}

if(FALSE) { # some tests
  bit_to_int(intToBits(-10))
  bit_to_int(intToBits(0))
  bit_to_int(intToBits(10))
}