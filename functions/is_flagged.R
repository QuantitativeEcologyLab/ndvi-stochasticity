#' to drop values flagged with quality flag assurance `.qa`. Note that the
#' flag positions are counted starting from 0 for consistency. For example,
#' if `flag = 2` and we want to exclude values (10, 11) and keep (00, 11),
#' we can use `QA %% (2 * 2) < 2`:
#' 00 keep; return TRUE  to keep
#' 01 keep; return TRUE  to keep
#' 10 drop; return FALSE to drop
#' 10 drop; return FALSE to drop
#' if you want to do the opposite, just negate the output with `!`
source('functions/bit_to_int.R')

is_flagged <- function(QA, flag_position) {
  if(length(QA) > 1) {
    stop('`QA` must be of length 1.\n')
  }
  if(length(flag_position) > 1) {
    stop('`flag_position` must be of length 1.\n')
  }
  bit <- numeric(16)
  bit[flag_position + 1] <- 1
  bit_value <- bit_to_int(bit = bit) 
  
  return(QA %% (bit_value * 2) >= bit_value)
}