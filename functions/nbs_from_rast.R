nbs_from_rast <- function(.r) {
  # get adjacent cells in all directions, and include the reference cell
  .nbs <-
    adjacent(.r, cells = cells(.r), directions = 8, include = TRUE) %>%
    as.data.frame() %>%
    transmute(
      ref_cell = V1, # first column is the starting cell
      # add the 8 surrounding neighbors
      adjacent = map(1:n(), \(i) {
        .z <- c(V2[i], V3[i], V4[i], V5[i], V6[i], V7[i], V8[i], V9[i])
        
        .values <- map_lgl(.z, \(.cell_id) {
          if(is.nan(.cell_id)) {
            return(NA)
          } else {
            return(.r[.cell_id]$z[1])
          }})
        
        .z <- .z[which(! is.na(.values))]
        
        if(length(.z) == 0) {
          return(0)
        } else {
          return(as.character(.z))
        }
      }))
  
  names(.nbs$adjacent) <- .nbs$ref_cell # add names of reference cells
  return(.nbs$adjacent) # return lists of reference cells only
}
