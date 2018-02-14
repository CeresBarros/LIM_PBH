## ------------------------------------------
## FUNCTIONS FOR FIRE_STSM MODULE
##
## Ceres: Dec 2017
## ------------------------------------------

## CROP & MASK TO STUDY AREA ----------------------------
## study.area and tocrop are "Raster*" or "Spatial*" objects
cropToStudyArea <- function(study.area, tocrop) {
  temp <- crop(x = tocrop, y = study.area)
  temp <- mask(x = temp, mask = study.area)
  return(temp)
}


## MAKE INITIAL VEGETATION RASTER IN STUDY AREA  --------
## area is a shapefile of the study area
## vegRas is a raster of vegetation types to be cropped amd reprojected
makePrefireVegetation <- function(area, vegRas){
  origCRS <- crs(area)
  area <- sp::spTransform(area, CRSobj = crs(vegRas))
  vegetation_prefire <- cropToStudyArea(study.area = area, tocrop = vegRas)
  
  vegetation_prefire <- projectRaster(from = vegetation_prefire, crs = origCRS, method = "ngb")
  
  return(vegetation_prefire)
}


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

## SEVERITY --------------------------------
## to be used with raster::calc
calculateSeverity <- function(x) {
  ## x is a vector of 3 values (pre-fire veg, post fire veg and fire presence)
  sev <- NA
  
  if(!all(is.na(x))) {
    ## throw error if vegetation changes without burning - this should not happen for now
    if(x[1] != x[2] & is.na(x[3])) {
      stop("Vegetation changes wihtout fire occurrence")
    }
    
    ## fire severity is > 0 if there was change
    if(x[1] != x[2]) {
      if(grepl("3|4|5", x[1])){
        ## highest severity when forests burn
        sev <- 3
      } else {
        stopifnot(x[1] == 2) ## check that pixel was initially a shrubland - otherwise something is wrong as hab0 and hab1 do not change
        ## if shrublands burn: medium severity
        sev <- 2
      }
    } else {
      ## low severity if grasslands burn
      if(x[1] == 1 & !is.na(x[3])) {
        sev <- 1
      } else {
        ## if there is no change, fire severity is 0
        sev <- 0
      }
    } 
  }
  
  return(sev)
}