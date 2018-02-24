## ------------------------------------------
## FUNCTIONS FOR simpleLCCSuccesion MODULE
##
## Ceres: Dec 2017
## ------------------------------------------

## VEGETATION TRANSITION --------------------

## Converts vegetation states, in burnt pixels, according to transition probability
## to be used with raster::calc
## fire_transitprobs is a transition matrix with n x n vegetation states, that needs to be in memory
## outputs a post-fire vegetation raster
vegTransition <- function(x) {
  temp <- x[1]
  if(!is.na(x[1]) & !is.na(x[2])) {
    temp <- names(which(fire_transitprobs[paste0("hab", temp), ] > runif(1)))
    ## convert to numeric
    temp <- as.numeric(sub("hab", "", temp))
  }
  return(temp)
}
