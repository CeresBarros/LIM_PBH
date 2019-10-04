## --------------------------------------
## Rsq FUNCTION FOR GAMLSSINF0TO1 MODELS
## --------------------------------------

## this is work around for the gamlss:::Rsq function which calculating Rsq in gamlssinf0to1 models
## because it doens't know which family to use.

Rsq_2 <- function (object, type = c("Cox Snell", "Cragg Uhler", "both")) {
  type <- match.arg(type)
  if (!is.gamlss(object))
    stop("this is design for gamlss objects only")
  Y <- if (object$family[1] %in% .gamlss.bi.list)
    cbind(object$y, object$bd - object$y) else object$y

  fam <- if (object$family == "InfBE")
    BEINF() else
    object$family

  suppressWarnings(m0 <- gamlssML(Y ~ 1, family = fam))
  rsq1 <- 1 - exp((2/object$N) * (logLik(m0)[1] - logLik(object)[1]))
  rsq2 <- rsq1/(1 - exp((2/object$N) * logLik(m0)[1]))
  if (type == "Cox Snell")
    return(rsq1)
  if (type == "Cragg Uhler")
    return(rsq2)
  if (type == "both")
    return(list(CoxSnell = rsq1, CraggUhler = rsq2))
}
