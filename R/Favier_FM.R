## Favier's model

rm(list = ls()); gc(reset = TRUE)
library(raster)
source("R/SpaDES/m/fire_spreadSTSM/R/FavierPercolationFM.R")

## modelled states: V, vegetated, F, burning, O, burnt and non-vegetated, E, empty
states <- factor(c("V", "F", "O", "E"), levels = c("V", "F", "O", "E"))

## parameter probabilities
# Q = probability of a vegetated cell igniting (V -> F)
# P = probability of extinction of a burning cell (F -> O)
# D = proportion of vegetated cells

D = 0.5 
Q = 0.5
P = 0.5

fireTransitions <- matrix(c(1-Q, Q, 0, 0,
                            0, 1-P, P, 0,
                            0, 0, 1, 0,
                            0, 0, 0, 1),
                          byrow = TRUE, nrow = 4, ncol = 4, 
                          dimnames = list(states, states))

## Favier's model - no markov chains
## lattice as raster
lattice <- raster(nrows= 102, ncols = 102); lattice[] <- NA
lattice[] <- apply(as.matrix(lattice[]), 1, FUN = function(x){
  if(runif(1) < D) return(which(states == "E")) else(return(which(states == "V")))  
})
lattice[1,] = lattice[102,] = lattice[,1] = lattice[,102] = NA

noFires = 1
steps = 50

startPix <- sample(which(!is.na(lattice[])), size = noFires, replace = FALSE)
lattice[startPix] <- which(states == "F")

neighbours <- data.table::data.table(adjacent(lattice, cells = which(!is.na(lattice[])), directions = 8, sorted = TRUE))

quickPlot::clearPlot()
quickPlot::dev()

lattice.ls <- list()
lattice.ls[[1]] <- lattice

for(k in 2:10){
  lattice.ls[[k]] <- FavierFM(lattice = lattice.ls[[k-1]])
  if(length(lattice.ls) == k){
    plot(lattice.ls[[k]], col = topo.colors(4), main = k)
  } else(break)
}
