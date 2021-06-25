#' Function to load raster stacks from .rds files that match a pattern
#'
#' @param x the pattern of file name to be matched
#' @param files the vector of all file names to be searched

loadStackFromRDS <- function(x, files) {
  stk <- lapply(grep(x, files, value = TRUE), readRDS) %>%
    stack(.)
  names(stk) <- renameFromFilenames(x, files)
  return(stk)
}

#' Function to load cohortData tables from .rds files that match a pattern
#'
#' @param x the pattern of file name to be matched
#' @param files the vector of all file names to be searched
#' @param pixelGroupMapStkList list of matching stacks of \code{pixelGroupMap}s
#'    List names must match \code{x}
#' @param yearSubset vector of years to subset

loadCohortDataFromRDS <- function(x, files, pixelGroupMapStk, yearSubset = NULL) {
  files <- grep(x, files, value = TRUE)

  pixelCohortData <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
    cohortData <- readRDS(ff)
    yr <- sub(".*year", "",  sub(".rds", "", ff))
    if (grepl("rep", ff)) {
      r <- sub(".*(rep)([0-9]+)/.*", "\\2", ff)
      rasName <- paste0("year", yr, "_rep", r)
    } else {
      r <- "1"
      rasName <- paste0("year", yr)
    }

    pixelGroupMap <- pixelGroupMapStk[[rasName]]
    cohortData <- addPixels2CohortData(cohortData, pixelGroupMap, doAssertion = FALSE) ## assertions where giving a very weird error about missing "sim"
    cohortData[, year := as.integer(yr)]
    cohortData[, rep := as.integer(r)]
    return(cohortData)
  }, pixelGroupMapStk = pixelGroupMapStk) %>%
    rbindlist(fill = TRUE, l = .)

  if (!is.null(yearSubset)) {
    pixelCohortData <- pixelCohortData[year %in% yearSubset]
  }

  return(pixelCohortData)
}


#' Function to make vegTypeData tables from vegTypeMap and pixelGroupMap stacks
#'
#' Assumes stack have "year" in the names. If they have "rep" that will be used too.
#'
#' @param vegTypeMapStk stack of \code{vegTypeMap}s
#' @param pixelGroupMapStk matching stacks of \code{pixelGroupMap}s
#'    List names must match \code{names(vegTypeMapStk)}
#' @param yearSubset vector of years to subset

vegTypeDataFromStks <- function(vegTypeMapStk, pixelGroupMapStk, yearSubset) {
  vegTypeSubset <- intersect(names(vegTypeMapStk), names(pixelGroupMapStk))

  vegTypeData <- lapply(vegTypeSubset, FUN = function(x) {
    yr <- grep("year", unlist(strsplit(x, split = "_")), value = TRUE)
    r <- if (grepl("rep", x)) {
      grep("rep", unlist(strsplit(x, split = "_")), value = TRUE)
    } else "1"

    data.table(pixelIndex = seq_len(ncell(vegTypeMapStk[[x]])),
               pixelGroup = getValues(pixelGroupMapStk[[x]]),
               vegType = vegTypeMapStk[[x]][],
               year = as.integer(sub("year", "", yr)),
               rep = as.integer(sub("rep", "", r)))
  }) %>%
    rbindlist(.)

  vegTypeData <- vegTypeData[!is.na(pixelGroup)]
  vegTypeData <- vegTypeData[year %in% yearSubset]
  return(vegTypeData)
}


#' Function to make pixelBurnData tables from rstCurrentFires stacks
#'
#' Assumes stack have "year" in the stack names. If they have "rep" that will be used too.
#'
#' @param rstCurrentFiresStk stack of \code{rstCurrentFires} maps

pixelBurnDataFromStks <- function(rstCurrentFiresStk) {
  pixelBurnData <- lapply(unstack(rstCurrentFiresStk), FUN = function(ras) {
    yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
    r <- if (any(grepl("rep", names(ras)))) {
      grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
    } else {
      1
    }
    data.table(pixelIndex = seq_len(ncell(ras)),
               burnt = as.integer(!is.na(getValues(ras))),
               fireID = as.integer(getValues(ras)),
               year = as.integer(sub("year", "", yr)),
               rep = as.integer(sub("rep", "", r)))
  }) %>%
    rbindlist(.)

  return(pixelBurnData)
}


#' Function to load severityData tables from .rds files that match a pattern
#'
#' Assumesfile names have "year". If they have "rep" that will be used too.
#'
#' @param x the pattern of file name to be matched
#' @param files the vector of all file names to be searched

loadSeverityDataFromRDS <- function(x, files) {
  files <- grep(x, files, value = TRUE)

  severityData <- lapply(files, FUN = function(ff) {
    severityData <- readRDS(ff)
    yr <- sub(".*year", "",  sub(".rds", "", ff))
    r <- if (grepl("rep", ff)) {
      sub(".*(rep)([0-9]+)/.*", "\\2", ff)
    } else "1"
    severityData[, year := as.integer(yr)]
    severityData[, rep := as.integer(r)]
    return(severityData)
  }) %>%
    rbindlist(fill = TRUE, l = .)

  return(severityData)
}

#' Function to get raster names from .rds files that match a pattern
#'
#' Assumes file names contain "year". If they contain "rep", this will be used too.
#' Assumes all files are .rds
#'
#' @param x the pattern of file names to be matched
#' @param files the vector of all file names to be searched

renameFromFilenames <- function(x, files) {
  if (any(grepl("rep", files))) {
    if (!all(grepl("rep", files))) {
      stop("Something is wrong. Some file names have 'rep', others don't")
    }
    paste(sub(".*year", "year", sub("\\.rds", "", grep(x, files, value = TRUE))),
          sub(".*(rep)([0-9]+)/.*", "\\1\\2",  grep(x, files, value = TRUE)), sep = "_")
  } else {
    sub(".*year", "year", sub("\\.rds", "", grep(x, files, value = TRUE)))
  }
}
