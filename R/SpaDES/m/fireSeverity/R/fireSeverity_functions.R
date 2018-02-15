## ------------------------------------------
## FUNCTIONS FOR fireSeverity MODULE
##
## Ceres: Dec 2017
## ------------------------------------------

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