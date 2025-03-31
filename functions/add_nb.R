library('dplyr')
library('sf')

add_nb <- function(p1, p2, add = FALSE) {
  if(! all(c(p1, p2) %in% names(nbs))) {
    stop('At least one of ', p1, ' or ', p2, ' is not in `names(nbs)`.')
  }
  
  # add neighbors if polygons are not already neighbors
  if(which(names(nbs) == p1) %in% nbs[[p2]] |
     which(names(nbs) == p2) %in% nbs[[p1]]) {
    stop('Polygons ', p1, ' and ', p2,
         ' are already neighbors! Check the values for each.')
  } else {
    filter(ecoregions, poly_id == p1 | poly_id == p2) %>%
      st_geometry() %>%
      st_transform('ESRI:102007') %>% # hawai'i albers
      plot(col = 1:2, main = paste('Polygons', p1, 'and', p2))
    
    if(add) {
      if(all(nbs[[p1]] == 0)) { # if no other neighbors do not add "0"
        nbs[[p1]] <<- which(names(nbs) == p2)
      } else { # otherwise keep other neighbors
        nbs[[p1]] <<- c(nbs[[p1]], which(names(nbs) == p2))
      }
      
      if(all(nbs[[p2]] == 0)) { # if no other neighbors do not add "0"
        nbs[[p2]] <<- which(names(nbs) == p1)
      } else { # otherwise keep other neighbors
        nbs[[p2]] <<- c(nbs[[p2]], which(names(nbs) == p1))
      }
      
      cat('Added', p1, 'as a neighbor to', p2, 'and vice-versa.\n')
    } else {
      warning('Did not add ', p1, ' as a neighbor to ', p2,
              ' and vice-versa!')
    }
  }
}
