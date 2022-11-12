#' Function to load raster stacks from .rds files that match a pattern
#'
#' @param x the pattern of file name to be matched
#' @param files the vector of all file names to be searched
#' @param startYear first year of the data, data files with an earlier year
#'   will be excluded
#' @param endYear last year of the data, data files with a later year
#'   will be excluded

loadStackFromRDS <- function(x, files, startYear, endYear) {
  filesYrs <- as.integer(sub(".*year", "",  sub(".rds", "", files)))
  files <- files[filesYrs >= startYear & filesYrs <= endYear]

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
#' @param startYear first year of the data, data files with an earlier year
#'   will be excluded
#' @param endYear last year of the data, data files with a later year
#'   will be excluded
#' @param yearSubset vector of years to subset

loadCohortDataFromRDS <- function(x, files, pixelGroupMapStk, startYear, endYear, yearSubset = NULL) {
  files <- grep(x, files, value = TRUE)

  filesYrs <- as.integer(sub(".*year", "",  sub(".rds", "", files)))
  files <- files[filesYrs >= startYear & filesYrs <= endYear]

  ## now subset to chosen years if need be
  if (!is.null(yearSubset)) {
    filesYrs <- as.integer(sub(".*year", "",  sub(".rds", "", files)))
    files <- files[filesYrs %in% yearSubset]
  }

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
#' @param startYear first year of the data, data files with an earlier year
#'   will be excluded
#' @param endYear last year of the data, data files with a later year
#'   will be excluded


loadSeverityDataFromRDS <- function(x, files, startYear, endYear) {
  files <- grep(x, files, value = TRUE)

  filesYrs <- as.integer(sub(".*year", "",  sub(".rds", "", files)))
  files <- files[filesYrs >= startYear & filesYrs <= endYear]

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

#' Wrapper to parallelize `convertToCNVegType`
#'
#' Assumes file names contain "year". If they contain "rep", this will be used too.
#' Assumes all files are .rds
#'
#' @param DT passed to `convertToCNVegType`
#' @param cPath used as `Cache(..., cacheRepo = cPath)`

parallelFUN <- function(DT, cPath) {
  rr <- .Random.seed
  a <- set.seed(123) # would be set.seed(P(sim)$.seed[["init"]]) )
  tempArg <- sample(1:nrow(DT), 100, replace = FALSE)
  .Random.seed <- rr # reset
  tempArg <- DT[tempArg,]
  setkey(DT, scenario, rep, year, pixelGroup)
  out <- Cache(convertToCNVegType,
               DT = DT,
               groupingCols = c("scenario", "rep", "year", "pixelGroup"),
               cachingArg = tempArg,
               omitArgs = c("DT"),
               cacheRepo = cPath,
               userTags = c("reportResults"))
  out
}

#' Classification of stand structure into Cameron Naficy's vegetation types.
#'
#' Uses a set of rules based on species relative biomass in a stand to
#' classify it  into one of 12 vegetation types:
#' "Oak", "PJ", "purePIPO", "DMCPIPO", "dryPSME", "PSME", "DMCPSME", "PICO",
#' "PIEN", "Broadleaf", "Mixedwood", "MMC".
#' This function deals with the eventuality of multiple matches between Cameron
#' species codes and simulated species codes
#'
#' @param DT A data.table with \code{groupingCols}, two columns with species codes
#'    (\code{speciesCode} and \code{Cameron}) and a column with realtive biomass per
#'    species (\code{relB}).
#' @param groupingCols character string of column names used for grouping and defining
#'    a single stand. Defaults to \code{"pixelGroup"}
#' @param pureCutoff threshold of relative biomass above which a stand is considered pure. Defaults to
#'    0.8.
#' @param drySp character spring of species characteristic of dry sites. Defaults to
#'    \code{c("PSME", "PIPO", "PIFL", "JUSC", "QUGA")}
#' @param moistSp character spring of species characteristic of moist sites. Defaults to
#'    \code{c("ABLA", "BEPA", "PIEN", "PIGL", "PIMO", "POBA", "THPL")}
#' @param cachingArg a vector of unique entries (e.g.) used to identify changes in \code{DT} when caching,
#'    used to avoid caching \code{DT}, when the object is large.
#'
#' @return a data.table with an extra column \code{"vegTypeCN"}
#' @export

convertToCNVegType <- function(DT, groupingCols = c("pixelGroup"), pureCutoff = 0.8,
                               drySp = c("PSME", "PIPO", "PIFL", "JUSC", "QUGA"),
                               moistSp = c("ABLA", "BEPA", "PIEN", "PIGL", "PIMO", "POBA", "THPL"),
                               cachingArg) {
  ## check:
  if (!all(c("speciesCode", "Cameron", groupingCols) %in% names(DT)))
    stop("not all groupingCols were found in DT") else
      setkeyv(DT, cols = groupingCols)   ## the function's not being applied by grouping.

  ## colums for the conversion
  cols <- c("Cameron", "speciesCode", "relB")

  ## Oak woodlands are dominated by oaks with no more dominant tree stature species
  DT[, oak := all(.sumRelBs("QUGA", .SD) >= pureCutoff,
                  .sumRelBs(c('PIPO', 'PSME', 'PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS'), .SD) < 0.05),
     by = groupingCols, .SDcols = cols]

  ## P-J woodlands are dominated by Pinyon juniper trees with no more dominant tree stature species
  DT[, PJ := all(.sumRelBs(c('PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS', 'QUGA'), .SD) >= pureCutoff,
                 .sumRelBs(c('PIPO', 'PSME'), .SD) < 0.05),
     by = groupingCols, .SDcols = cols]

  ## Pure PIPO if PIPO is heavily dominant and accompanied by small amount of other species
  DT[, purePIPO := all(.sumRelBs("PIPO", .SD) >= pureCutoff,
                       .sumRelBs(c('PSME', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'), .SD) < 0.30),
     by = groupingCols, .SDcols = cols]

  ## DMC if PIPO present at >= 10% but less than 70% and other species are all dry site species
  DT[, DMCPIPO := all(.sumRelBs("PIPO", .SD) >= 0.10,
                      .sumRelBs("PIPO", .SD) < pureCutoff,
                      .sumRelBs(c('PSME', 'PIPO', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'), .SD) >= 0.50),
     by = groupingCols, .SDcols = cols]


  ## If PSME is dominant and dry site species are present
  DT[, dryPSME := all(.sumRelBs("PSME", .SD) >= pureCutoff,
                      .sumRelBs(moistSp, .SD) < 0.10,
                      .sumRelBs(c('JUSC', 'JUOC', 'JUOS', 'PIFL', 'PIED', 'PIMO2', 'QUGA'), .SD) > 0.05),
     by = groupingCols, .SDcols = cols]

  ## If PSME is dominant and dry site species are absent
  DT[, PSME := all(.sumRelBs("PSME", .SD) >= pureCutoff,
                   .sumRelBs(moistSp, .SD) < 0.10,
                   .sumRelBs(c('JUSC', 'JUOC', 'JUOS', 'PIFL', 'PIED', 'PIMO2', 'QUGA'), .SD) <= 0.05),
     by = groupingCols, .SDcols = cols]


  ## DMC if ponderosa pine not present, dry site species are dominant but may be
  ## micov.BAed with some other species (e.g. POTR, LAOC, PICO), and few moist site species are present in significant numbers
  DT[, DMCPSME := all(.sumRelBs(drySp, .SD) >= 0.50,
                      .sumRelBs(moistSp, .SD) < 0.10),
     by = groupingCols, .SDcols = cols]

  ## PICO if PICO dominates stand
  DT[, PICO := .sumRelBs("PICO", .SD) >= 0.5,
     by = groupingCols, .SDcols = cols]

  ## lowland PIEN if dominated by spruce
  DT[, PIEN := .sumRelBs(c('PIEN', 'PIEN/PIGL', 'PIGL', 'ABLA'), .SD) > 0.50,
     by = groupingCols, .SDcols = cols]

  ## broadleaf and mixedwood
  DT[, broadleaf := .sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), .SD) >= pureCutoff,
     by = groupingCols, .SDcols = cols]
  DT[, mixedwood := all(.sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), .SD) < pureCutoff,
                        .sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), .SD) >= 0.25),
     by = groupingCols, .SDcols = cols]

  ## MMC if all others are false
  cols <- setdiff(names(DT), c("speciesCode", "Cameron", "relB", groupingCols))
  DT[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
  DT[, MMC := rowSums(.SD) == 0, .SDcols = cols]

  toResolve <- DT[MMC != TRUE]
  resolved <- DT[MMC == TRUE]

  ## colums must be in this order to resolve the veg type
  ## the first to be TRUE is the final veg type (MMC is ignored)
  colsVeg <- c("oak", "PJ", "purePIPO", "DMCPIPO", "dryPSME", "PSME", "DMCPSME", "PICO", "PIEN", "broadleaf", "mixedwood")
  colsNotVeg <- setdiff(names(toResolve), colsVeg)
  allCols <- c(colsNotVeg, colsVeg)
  toResolve <- toResolve[, ..allCols]  ## reorder.

  toResolve[, vegTypeCN := apply(as.matrix(.SD), 1, which.max), .SDcols = colsVeg]
  toResolve[, vegTypeCN := colsVeg[vegTypeCN]] ## replace by col name

  ## now convert the MMC col to the vegTypeCN
  resolved[, vegTypeCN := "MMC"]

  ## subset and join
  colsNotVeg <- setdiff(colsNotVeg, "MMC")
  allCols <- c(colsNotVeg, "vegTypeCN")

  rbind(toResolve[, ..allCols], resolved[, ..allCols], use.names = TRUE)
}

#' internal function that sums relative biomasses for species matching a character string,
#' but that can be appear duplicated in another species coding column.
#' @param sppToMatch character string of species to match against for summing B.
#' @param DT data.table with columns 'Cameron', 'speciesCode', 'relB'.
.sumRelBs <- function(sppToMatch, DT) {
  DT[Cameron %in% sppToMatch] %>%
    .[, .(speciesCode, relB)] %>%
    unique(.) %>%
    .$relB %>%
    sum(.)
}
