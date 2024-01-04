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

pixelBurnDataFromStks <- function(rstCurrentFiresStk, i = NULL) {
  if (is.null(i)) i <- names(rstCurrentFiresStk)

  pixelBurnData <- lapply(unstack(rstCurrentFiresStk[[i]]), FUN = function(ras) {
    yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
    r <- if (any(grepl("rep", names(ras)))) {
      grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
    } else {
      1
    }

    if (inherits(ras, "Raster")) ras <- rast(ras)

    DT <- as.data.table(ras, cell = TRUE)
    setnames(DT, c("pixelIndex", "fireID"))
    DT[, burnt := as.integer(!is.na(fireID))]
    DT[, year := as.integer(sub("year", "", yr))]
    DT[, rep := as.integer(sub("rep", "", r))]
  }) %>%
    rbindlist(.)

  return(pixelBurnData)
}

#' Function to make a data.table from yearly stacks or rasters
#' of fire attributes
#'
#' Assumes list names are <year>.<rep>, or just <year>
#'
#' @param fireAttrRasLs list of stacks of fire attribute maps
#'
#' @return data.table of raster attributes with columns
#'  "year", "rep" and "pixelIndex"
#'
#' @importFrom data.table rbindlist

fireAttrDTFromRasLs <- function(fireAttrRasLs, i = NULL) {
  if (is.null(i)) i <- names(fireAttrRasLs)
  DTLs <- Map(i = i,
              f = .getFireAttr,
              MoreArgs = list(fireAttrRasLs = fireAttrRasLs))
  rbindlist(DTLs, use.names = TRUE)
}

#' @importFrom data.table as.data.table setnames
.getFireAttr <- function(i, fireAttrRasLs) {
  ii <- sub("year", "", sub("_rep", ".", i, fixed = TRUE))
  yr <- strsplit(ii, split = "\\.")[[1]][1]
  r <- strsplit(ii, split = "\\.")[[1]][2]
  if (is.na(r)) {
    r <- 1
  }

  DT <- as.data.table(fireAttrRasLs[[i]], cell = TRUE)
  setnames(DT, "cell", "pixelIndex")
  DT[, `:=`(year = as.integer(yr),
            rep = as.integer(r))]
  DT
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

#' Makes a stack raster of severity class and actual severity (B lost)
#'
#' @param rasterToMatch raster layer to use to rasterize fire data.
#'   Must have `ncell` matching `pixelIndex` in `sevData`
#' @param sevData data.table with columns listed in `sevCols`.
#' @param sevCols character. Names of columns with severity data to extract to
#'   make severity rasters. Column names will determine raster names.
#'
#' @return SpatRaster stack of raster of severity properties defined
#'  in `sevCols`.
#'
#' @export
makeSevRasters <- function(rasterToMatch, sevData,
                           sevCols = c("severity", "severityB", "severityPropB")) {
  if (!"severity" %in% names(sevData)) {
    sevData[, severity := 5L]
  }

  RTM <- rasterToMatch
  if (inherits(RTM, "Raster")) {
    RTM <- rast(RTM)
  }

  RTM[] <- NA_integer_

  sevRasters <- sapply(sevCols, function(ccol, sevData, RTM) {
    sevRaster <- RTM

    sevRaster[sevData[["pixelIndex"]]] <- sevData[[ccol]]

    ## for now add a 0 severity (class and B) in unburnt pixels
    RTM[!is.na(as.vector(RTM[]))] <- 0L
    sevRaster <- terra::cover(sevRaster, RTM, values = NA)

    names(sevRaster) <- ccol
    sevRaster
  }, sevData = sevData, RTM = RTM, simplify = FALSE, USE.NAMES = TRUE)

  return(rast(sevRasters))
}

#' Calculates patch size as the area in log-ha of discrete
#' patches of distinct severity classes. All fire patches
#' (in and outside of forested pixels) are considered. Patches
#' in non-forested pixels are assumed to have 1 severity class
#' (`severity`) with 0 lost biomass (`severityB` and `severityBProp`).
#'
#' Adapted from [Steel et al 2021](https://github.com/zacksteel/pyrodiversity/blob/main/code/patch_surface.R).
#'
#' @param sevClassRas a raster of severity (in class).
#' @export
calcPatchSize <- function(sevClassRas, fireRas) {
  ## check that all pixels with severity are in a fire perimeter
  ## the inverse does not need to be true (non-forested pixels can burn, but may not have severity)
  pixWSev <- which(as.vector(sevClassRas[]) > 0)
  if (any(is.na(fireRas[pixWSev]))) stop("Found pixels with severity outside of a fire perimeter")

  ## convert to terra if need be
  if (inherits(fireRas, "Raster")) {
    fireRas <- rast(fireRas)
  }

  if (inherits(sevClassRas, "Raster")) {
    sevClassRas <- rast(sevClassRas)
  }

  ## loop through fire IDs
  fireIDs <- na.exclude(unique(as.vector(fireRas[])))

  ## storage raster
  patchSizeRas <- rast(sevClassRas)

  if (length(fireIDs)) {
    for (fire in fireIDs) {
      focalFireRas <- fireRas
      focalFireRas[as.vector(focalFireRas[]) != fire] <- NA

      ## mask severity to fire perimeter -- pixels outside fire perimeter will get an NA for patch size.
      focalFireSevRas <- mask(sevClassRas, focalFireRas)

      ## convert to vector
      focalFireSevPoly <- as.polygons(focalFireSevRas, dissolve = TRUE)
      focalFireSevPoly <- disagg(focalFireSevPoly)

      ## Calculate log-area
      focalFireSevPoly$log_area <- log(expanse(focalFireSevPoly, unit = 'ha'))

      ## rasterize and add to output raster
      focalFireSevRas <- rasterize(focalFireSevPoly, patchSizeRas, field = "log_area")

      patchSizeRas <- cover(patchSizeRas, focalFireSevRas)
    }
  } else {
    patchSizeRas[] <- NA_integer_
    names(patchSizeRas) <- "log_area"
  }

  patchSizeRas
}


#' Assesses with pixels have no fire history from a stack
#' of yearly fire perimeter rasters
#'
#' @param rstCurrentFiresStk stack of named yearly rasters of fire perimeters.
#'   If `r` is a numeric and supplied, the function will look for "rep<r>"
#'   in the names to subset the raster layers.
#' @param r integer. The replicate number to subset raster layers to.
#'   If `NULL` the function returns a `data.table` with column `rep == 1L`
#' @param rasterToMatch raster delimiting pixels used for simulation.
#'
#' @return `data.table` of unburnt pixels with columns `pixelIndex`, `burnt`
#'  (only with 0L), and `rep`
#' @export
makeNoFireHistoryData <- function(rstCurrentFiresStk, rasterToMatch, r = NULL,
                                  doAssertion = getOption("LandR.assertions", TRUE)) {
  ## covert to terra
  if (inherits(rstCurrentFiresStk, "Raster")) {
    rstCurrentFiresStk <- rast(rstCurrentFiresStk)
  }

  if (inherits(rasterToMatch, "Raster")) {
    rasterToMatch <- rast(rasterToMatch)
  }


  if (!is.null(r)) {
    ## subset to rep
    subsetRas <- grep(paste0("rep", r), names(rstCurrentFiresStk))
    rstCurrentFiresStk <- rstCurrentFiresStk[[subsetRas]]
  } else {
    r <- 1L
  }

  ## make a raster with 0/1 fore pixels w/o or w/ fire history, respectively
  burnRas <- sum(rstCurrentFiresStk, na.rm = TRUE)
  burnRas[as.vector(is.na(burnRas[]))]  <- 0L  ## NAs have no fire history

  ## now mask again to put NAs outside SA
  burnRas <- mask(burnRas, rasterToMatch)

  ## pixels with > 0 have a fire history
  burnRas[as.vector(burnRas[]) > 0] <- 1L
  names(burnRas) <- "burnt"

  ## calculate unburnt patch area - as in Steel et al 2021.
  burnRasPoly <- as.polygons(burnRas, dissolve = TRUE)
  burnRasPoly <- disagg(burnRasPoly)

  ## Calculate log-area
  burnRasPoly$log_area <- log(expanse(burnRasPoly, unit = 'ha'))

  ## rasterize and add to output raster
  patchSizeRas <- rasterize(burnRasPoly, burnRas, field = "log_area")

  ## calculate patch size in pixels
  burnRasPoly$ID <- 1:length(burnRasPoly)

  unburnPatches <- rasterize(burnRasPoly, burnRas, field = "ID")
  patchSizeNoPix <- as.data.table(unburnPatches, cell = TRUE)
  patchSizeNoPix[, patchSizePix := .N, by = ID]

  ## make final DT with patch sizes
  noFireHistoryData <- as.data.table(patchSizeRas, cell = TRUE)
  noFireHistoryData <- noFireHistoryData[patchSizeNoPix, on = "cell"]
  setnames(noFireHistoryData, c("cell", "log_area"), c("pixelIndex", "patchSizeLogHa"))

  ## check
  if (doAssertion) {
    test <- setdiff(which(!is.na(as.vector(rasterToMatch[]))), noFireHistoryData$pixelIndex)
    if (length(test)) {
      stop("All active pixels in study area should exist in this table; please debug 'makeNoFireHistoryData'.")
    }
  }

  unburntPix <- which(as.vector(burnRas[]) == 0)
  noFireHistoryData <- noFireHistoryData[pixelIndex %in% unburntPix]
  noFireHistoryData[, rep := r]

  return(noFireHistoryData)
}

#### OLD EVENT:
# loadSimulationDataEvent <- function(sim) {
#   # ! ----- EDIT BELOW ----- ! #
#   mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...
#
#   grepPattrn <- paste0("/", P(sim)$scenarios)
#   if (all(!is.na(P(sim)$reps))) {
#     reps <- P(sim)$reps
#     reps[reps < 10] <- paste0("0?", reps[reps < 10])
#     grepPattrn <- paste0(grepPattrn, "_",
#                          "rep", "(", paste0(reps, collapse = "|"), ")/")
#   }
#   names(grepPattrn) <- P(sim)$scenarios
#
#   ## make list of stacked rasters -- filters between start and end years
#   rstCurrentFiresStkList <- sapply(grepPattrn,
#                                    FUN = loadStackFromRDS,
#                                    files = sim$rstCurrentFiresFiles,
#                                    startYear = P(sim)$startYear,
#                                    endYear = P(sim)$endYear,
#                                    simplify = FALSE, USE.NAMES = TRUE)
#
#   pixelGroupMapStkList <- sapply(grepPattrn,
#                                  FUN = loadStackFromRDS,
#                                  files = sim$pixelGroupMapFiles,
#                                  startYear = P(sim)$startYear,
#                                  endYear = P(sim)$endYear,
#                                  simplify = FALSE, USE.NAMES = TRUE)
#
#   vegTypeMapStkList <- sapply(grepPattrn,
#                               FUN = loadStackFromRDS,
#                               files = sim$vegTypeMapFiles,
#                               startYear = P(sim)$startYear,
#                               endYear = P(sim)$endYear,
#                               simplify = FALSE, USE.NAMES = TRUE)
#
#
#   ## pixelCohortData tables -- filters between start and end years and yearSubset
#   pixelCohortDataList <- mapply(FUN = loadCohortDataFromRDS,
#                                 x = grepPattrn,
#                                 pixelGroupMapStk = pixelGroupMapStkList,
#                                 MoreArgs = list(files = sim$cohortDataFiles,
#                                                 yearSubset = P(sim)$yearSubset,
#                                                 startYear = P(sim)$startYear,
#                                                 endYear = P(sim)$endYear),
#                                 SIMPLIFY = FALSE, USE.NAMES = TRUE)
#
#   ## vegTypeData tables -- filters to yearSubset
#   vegTypeDataList <- mapply(FUN = vegTypeDataFromStks,
#                             vegTypeMapStk = vegTypeMapStkList,
#                             pixelGroupMapStk = pixelGroupMapStkList,
#                             MoreArgs = list(yearSubset = P(sim)$yearSubset),
#                             SIMPLIFY = FALSE, USE.NAMES = TRUE)
#
#   ## pixelBurnData tables - all rasters
#   pixelBurnDataList <- sapply(rstCurrentFiresStkList,
#                               FUN = pixelBurnDataFromStks,
#                               simplify = FALSE, USE.NAMES = TRUE)
#
#   ## add scenario column when binding
#   ## exclude NAs early to save memory when binding
#   pixelBurnDataList <- lapply(pixelBurnDataList, function(allPixelBurnData) allPixelBurnData[fireID != "NA"])
#   allPixelBurnData <- rbindlist(pixelBurnDataList, idcol = "scenario", use.names = TRUE)
#
#   ## severityData tables -- filters between start and end years
#   severityDataList <- sapply(grepPattrn,
#                              FUN = loadSeverityDataFromRDS,
#                              files = sim$severityDataFiles,
#                              startYear = P(sim)$startYear,
#                              endYear = P(sim)$endYear,
#                              simplify = FALSE, USE.NAMES = TRUE)
#   ## add scenario column when binding
#   allSeverityData <- rbindlist(severityDataList, idcol = "scenario", fill = TRUE, use.names = TRUE)
#
#   ## join fire data
#   amc::.gc()
#   allPixelBurnData <- allSeverityData[allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]
#
#   ## clean ws
#   rm(rstCurrentFiresStkList, vegTypeMapStkList, pixelGroupMapStkList,
#      pixelBurnDataList, severityDataList, allSeverityData)
#   amc::.gc()
#
#   ## join tables, add scenario col and rbind.
#   ## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
#   ## pixelBurntData can also have a different number of years if the saving frequency differs
#   ## so join by keeping all pixels, calculate fire properties per pixel, then subset to
#   ## cohort data years.
#   pixelCohortDataList <- mapply(FUN = merge,
#                                 x = vegTypeDataList,
#                                 y = pixelCohortDataList,
#                                 MoreArgs = list(by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE),
#                                 SIMPLIFY = FALSE, USE.NAMES = TRUE)
#   amc::.gc()
#
#   ## add scenario column when binding
#   allPixelCohortData <- rbindlist(pixelCohortDataList, use.names = TRUE, idcol = "scenario")
#   if (exists("allPixelCohortData")) rm(pixelCohortDataList, vegTypeDataList)
#   amc::.gc()
#
#   ## checks
#   if (mod$doAssertion)  {
#     ## note that some years may not have all the reps in allPixelBurnData if there where no fires
#     ## because we excluded NAs (i.e. pixels without fires)
#     repsFire <- unique(allPixelBurnData[, length(unique(rep)), by = .(scenario)]$V1)
#     repsCohortData <- unique(allPixelCohortData[, length(unique(rep)), by = .(scenario, year)]$V1)
#     test1 <- all(length(repsFire) == 1, length(repsCohortData) == 1)
#     test2 <- all(length(repsFire) == length(repsCohortData),
#                  identical(sort(repsFire), sort(repsCohortData)))
#
#     if (!isTRUE(test1)) {
#       stop("Fire and/or cohort data do not have the same number of reps across scenarios")
#     }
#     if (!isTRUE(test2)) {
#       stop("Fire and cohort data differ in number of reps (per scenario/year)")
#     }
#   }
#
#   ## export to sim
#   sim$allPixelBurnData <- allPixelBurnData
#   sim$allPixelCohortData <- allPixelCohortData
#
#   # ! ----- STOP EDITING ----- ! #
#   return(invisible(sim))
# }

