#' Calculate fire weather indices from weather table
#'
#' This is a wrapper function ued to cache fire weather indices
#'   calculations from a weather table. See \code{cffdr::fwi}
#'   for further information.
#'
#' @param tableSize A vector of row IDs used for caching purposes only
#'    to avoid caching the \code{fwi} \code{input} argument, which can
#'    be a very large table
#' @param ... arguments passed to \code{fwi}
#' @param suppressWarnings Defaults to \code{TRUE} to suppress warnings
#'    from \code{fwi}.
#'
#' @return a table with fire weather indices. See \code{cffdr::fwi} for
#'    further information
#'
#' @seealso \code{cffdr::fwi}
#' @export
#' @importFrom cffdrs fwi

calcFWIWeatherTable <- function (tableSize, ...,
                                 suppressWarnings = TRUE) {
  FWIoutputs <- if (suppressWarnings) {
    suppressWarnings({
      fwi(...)
    })
  } else fwi(...)

  return(FWIoutputs)
}
