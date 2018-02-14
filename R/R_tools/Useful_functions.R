## ------------------------------------------------------
## USEFUL FUNCTIONS
##
## Ceres: Dec 2017
## ------------------------------------------------------

## this script should be sourced

## CHECK PROJECTIONS ------------------------------------
## obj.list is a list of spatial objects
checkProjections <- function(obj.list){
  projs <- sapply(obj.list, FUN = function(x) {
    return(
      eval(expr = parse(text = paste0("projection(", x, ")")))
    )
  })
  return(projs)
}


## CROP & MASK TO STUDY AREA ----------------------------
## study.area and tocrop are "Raster*" or "Spatial*" objects
cropToStudyArea <- function(study.area, tocrop) {
  temp <- crop(x = tocrop, y = study.area)
  temp <- mask(x = temp, mask = study.area)
  return(temp)
}


## RASTER TO BINARY MATRIX ------------------------------
## converts a vector of values into a binary "presence/absence" matrix
## x is a vector
vector2binmatrix <- function(x) {
  x <- as.character(x)
  return(model.matrix( ~ x-1))
}


## FIND NEIGHBOURS IN MATRIX ----------------------------
## finds the 8 neighbours of each cell in a matrix
## output is a matrix of 8 rows, with a columns per cell of the input matrix, 
## which is treated by columns

neighboursMatrix = function(mat) {
  mat2 <- cbind(NA, rbind(NA, mat ,NA), NA)   ## makes a border of NAs
  addresses <- expand.grid(x = 1:nrow(mat), y = 1:ncol(mat)) ## all matrix coordinates
  neighs <- c()
  for(i in 1:-1) {
    for(j in 1:-1) {
      if(i != 0 || j != 0) {
        neighs <- rbind(neighs,mat2[addresses$x+i+1+nrow(mat2)*(addresses$y+j)])   ## each column contains the neighboors a cell (going by columns in mat) 
      }
    }
  }
  return(neighs)
}
