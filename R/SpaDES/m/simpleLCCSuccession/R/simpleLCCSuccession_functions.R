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
  ## if both have values, then its a burnt veg. pixel
  ## climate (x[3]) increases the probability of burning the vegetation
  if(!is.na(x[1]) & !is.na(x[2])) {
    temp <- names(which(fireTransitMatrix[paste0("hab", temp), ] > runif(1)))
    ## convert to numeric
    temp <- as.numeric(sub("hab", "", temp))
  }
  
  return(temp)
}

vegTransition2 <- function(x) {
  temp <- x[1]
  
  ## if both have values, then its a burnt veg. pixel
  ## climate (x[3]) increases the probability of burning the vegetation
  if(!is.na(x[1]) & !is.na(x[2])) {
    temp <- names(which(fireTransitMatrix[paste0("hab", temp), ] * x[3] > runif(1)))
    if(length(temp) > 0) {
      ## convert to numeric
      temp <- as.numeric(sub("hab", "", temp))
    } else {
      temp <- x[1]
    }
  } else if(!is.na(x[1])) {
    ## no fire - vegetation succession procedes
      temp <- names(which(vegTransitMatrix[paste0("hab", temp), ] > runif(1)))
      ## convert to numeric
      temp <- as.numeric(sub("hab", "", temp))
    }
  
  return(temp)
}
