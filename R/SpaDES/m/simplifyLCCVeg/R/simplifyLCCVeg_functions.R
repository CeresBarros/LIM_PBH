## ------------------------------------------
## FUNCTIONS FOR simplifyLCCVeg MODULE
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