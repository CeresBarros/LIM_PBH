## CLEANING DAVE's FIRE DATA

library(raster)

AL1_prefire <- shapefile("data/fires_Dave/Projected_renamed/albertafires1_prefire")
AL2_prefire <- shapefile("data/fires_Dave/Projected_renamed/albertafires2_prefire")
AL1_postfire <- shapefile("data/fires_Dave/Projected_renamed/albertafires1_postfire")
AL2_postfire <- shapefile("data/fires_Dave/Projected_renamed/albertafires2_postfire")

## which variables only have NAs?
varNAs <- names(which(sapply(lapply(AL1_prefire@data, FUN = function(x) unique(x)),
       FUN = function(x) all(is.na(x)))))

## which variables only have 0s?
var0s <- names(which(sapply(lapply(AL1_prefire@data[AL1_prefire@data[, !names(AL1_prefire@data) %in% varNAs]], 
                                   FUN = function(x) {
                                     if(class(x) != "character") sum(x)}),
                            FUN = function(x) any(x == 0))))
## show unique values of character variables
lapply(AL1_prefire@data, FUN = function(x) {
  if(class(x) == "character") unique(x)
})

## alberta shapefiles have many seemingly duplicated variables
dups <- sub("_$", "", names(AL1_prefire@data))
dups <- dups[duplicated(dups)]

## of these are different?
sapply(dups, FUN = function(j){
  temp <- AL1_prefire@data[, which(sub("_$", "", names(AL1_prefire@data)) == j)]
  # temp[is.na(temp)] <- "-99999"
  any(temp[,1] != temp[,2], na.rm = TRUE)
})
## a priori, variables with "_" are more complete. To check with Dave.    

all(sub(" .*", "", AL2_prefire@data$FORESTKEY) == AL2_prefire@data$FIRE_NAME, na.rm = TRUE)


##


