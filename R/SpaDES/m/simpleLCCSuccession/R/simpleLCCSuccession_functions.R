## ------------------------------------------
## FUNCTIONS FOR simpleLCCSuccesion MODULE
##
## Ceres: Dec 2017
## ------------------------------------------

## VEGETATION TRANSITION --------------------

## Converts vegetation states, in burnt pixels, according to a transition probability
## to be used with raster::calc
## x is a two raster stack with the vegetation and fire rasters.
## transitionMatrix is a transition matrix with n by n vegetation states - needs to be in memory, 
##as calc doesn't not pass extra arguments to functions

## outputs a post-fire vegetation raster

vegTransition <- function(x) {
  temp <- x[1]
  if(!is.na(x[1]) & !is.na(x[2])) {
    temp <- names(which(transitionMatrix[paste0("hab", temp), ] > runif(1)))
    ## convert to numeric
    temp <- as.numeric(sub("hab", "", temp))
  }
  return(temp)
}
