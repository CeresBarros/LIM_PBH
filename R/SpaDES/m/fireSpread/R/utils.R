## ------------------------------------------------------
## USEFUL FUNCTIONS
##
## Ceres: Sep 2018
## ------------------------------------------------------

## FUNCTION TO CALCULATE THE MODE
## x is a vector of values, numeric or character/factor
## ... can be used to pass arguments to table and max
getmode <- function(x, ...) {
  tab <- table(x)
  if (is.numeric(x)) {
    as.numeric(names(which(tab == max(tab))))
  } else {
    as.character(names(which(tab == max(tab))))
  }
}
