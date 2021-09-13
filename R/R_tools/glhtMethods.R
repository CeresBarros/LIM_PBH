## GLS METHODS FOR glht

model.matrix.gls <- function(object, ...)
  model.matrix(terms(object), data = getData(object), ...)


model.frame.gls <- function(object, ...)
  model.frame(formula(object), data = getData(object), ...)


terms.gls <- function(object, ...)
  terms(model.frame(object),...)
