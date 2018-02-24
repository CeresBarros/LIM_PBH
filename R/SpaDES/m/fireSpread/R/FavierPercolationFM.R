## ------------------------------------------
## FAVIERS PERCOLATION FIRE MODEL
##
## Ceres: Feb 2018
## ------------------------------------------

## see origianl paper by Favier C (2004) Percolation model of fire dynamic. Physics Letters, Section A: General, Atomic and Solid State Physics, 330, 396–401.

FavierFM <- function(lattice){
  ## get fire cell positions
  fireCells <- which(lattice[] == which(states == "F"))
  
  if(length(fireCells) > 0){
    ## focus on neighbours of a fire cell - get their positions
    fireNeighbours <- data.table::data.table(adjacent(lattice, 
                                                      cells = which(lattice[] == which(states == "F")), 
                                                      directions = 8, sorted = TRUE))
    fireNeighbours <- setdiff(fireNeighbours$to, which(is.na(lattice[]))) ## remove NAs if any
    
    
    ## spread fire to adjacent cells
    temp <- apply(as.matrix(fireNeighbours), 1,  FUN = function(i){
      x <- lattice[][i]   ## get cell value
      xx = NA
      
      ## get cell values surrounding each neighbour
      surroundingVals <- lattice[adjacent(lattice, cells = i, directions = 8, sorted = TRUE)[,2]]
      
      ## Keep focal cell vegetated, or set it on fire according to
      ## calculate the probability of fire spreading to focal cell
      if(x == which(states == "V")){
        xx <- which(states == "V")
        pBurning <- 1-(1-Q)^sum(surroundingVals == which(states == "F"), na.rm = TRUE)
        if(runif(1) < pBurning) xx <- which(states == "F")
      }
      ## Burnt (O) and empty (E) cells remain as such
      if(x == which(states == "O")) xx = x
      if(x == which(states == "E")) xx = x
      return(xx)
    })
    
    ## extinguish fires started before
    temp2 <- apply(as.matrix(fireCells), 1, FUN = function(i){
      x = lattice[][i]
      xx = NA
      if(x == which(states == "F")){
        if(runif(1) < P) xx <- which(states == "O") else xx <- x
      } else(stop("fireCell is not of fire... something's wrong"))
      return(xx)
    })
    
    lattice[fireNeighbours] <- temp
    lattice[fireCells] <- temp2
    
    return(lattice)
  } 
}