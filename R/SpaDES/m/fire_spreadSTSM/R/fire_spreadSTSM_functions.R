## ------------------------------------------
## FUNCTIONS FOR fire_spreadSTSM MODULE
##
## Ceres: Dec 2017
## ------------------------------------------

## CALCULATE FIRE SPREAD RASTER ------------------------
## pixels is the number of pixels where fires will start 
## fMask is a raster of burnable areas (mask format)
## fROS is a raster of Rats Of Spread 
## fSize is the final fire size in no. pixels
makeFireSpreadRas <- function(pixels, fMask, fROS, fSize){
  ## generate starting fire pixels
  startPix <- sample(which(!is.na(fMask[])), pixels)
  
  ## calculate spread
  spreadRas <- spread2(landscape = fMask, spreadProbRel = fROS, 
                       start = startPix, maxSize = fSize, plot.it = FALSE)
  
  ## remove fires outside burnable areas
  spreadRas[is.na(fMask)] <- NA
  
  ## change it to mask format
  # spreadRas[!is.na(spreadRas)] <- 1
  return(spreadRas)
}

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

