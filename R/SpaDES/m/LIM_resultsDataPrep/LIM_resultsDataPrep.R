## Everything in this file and any files in the R directory are sourced during `simInit()`;
## all functions and objects are put into the `simList`.
## To use objects, use `sim$xxx` (they are globally available to all modules).
## Functions can be used inside any function that was sourced in this module;
## they are namespaced to the module, just like functions in R packages.
## If exact location is required, functions will be: `sim$.mods$<moduleName>$FunctionName`.
defineModule(sim, list(
  name = "LIM_resultsDataPrep",
  description = "Module to compile results from LIM model simulations",
  keywords = c("results", "LIM project", "data prep."),
  authors = structure(list(list(given = "Ceres", family = "Barros", role = c("aut", "cre"),
                                email = "cbarros@mail.ubc.ca", comment = NULL)), class = "person"),
  childModules = character(0),
  version = list(SpaDES.core = "1.1.0.9004",
                 LIM_resultsDataPrep = "0.0.1"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.md", "LIM_resultsDataPrep.Rmd")), ## same file
  reqdPkgs = list("crayon", "data.table", "dplyr", "future", "future.apply",
                  "raster", "terra",
                  "PredictiveEcology/LandR@development (>= 1.0.7.9023)",
                  "PredictiveEcology/reproducible@development (>= 1.2.11)",
                  "PredictiveEcology/SpaDES.core@development (>= 1.1.0.9004)",
                  "CeresBarros/ToolsCB"),
  parameters = rbind(
    defineParameter("endYear", "integer", 2111L, NA_integer_, NA_integer_,
                    "The last year of simulation results to use."),
    defineParameter("ncores", "integer", 8L, 1L, NA_integer_,
                    "Number of cores to use if P(sim)$parallel is TRUE"),
    defineParameter("parallel", "logical", TRUE, NA, NA,
                    paste("Should data processing be parallelized? Currently only used in assigning",
                          "functional vegetation types following Cameron Naficy's classification")),
    defineParameter("reps", "integer", 1L:10L, NA_integer_, NA_integer_,
                    "The simulation repetitions to compile. If no repetitions were performed set to NA"),
    defineParameter("scenarios", "character", c("PM", "noPM"), NA, NA,
                    paste("The simulation scenarios to compile - must correspond with the names used in the",
                          "outputs folder tree")),
    defineParameter("startYear", "integer", 2011L, NA_integer_, NA_integer_,
                    "The first year simulation results to use"),
    defineParameter("yearSubset", "integer", as.integer(unique(c(seq(2011, 2111, 5), 2111))), NA_integer_, NA_integer_,
                    paste("Specific simulation years to compile - only vegetation dynamics will be subset.",
                          "outputs folder tree. If using all years set to NULL. Must contain `startYear` and `endYear`")),
    defineParameter(".plots", "character", "screen", NA, NA,
                    "Used by Plots function, which can be optionally used here"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "Describes the simulation time at which the first plot event should occur."),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "Describes the simulation time interval between plot events."),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "Describes the simulation time at which the first save event should occur."),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events."),
    ## .seed is optional: `list('init' = 123)` will `set.seed(123)` for the `init` event only.
    defineParameter(".seed", "list", list(), NA, NA,
                    "Named list of seeds to use for each event (names)."),
    defineParameter(".useCache", "logical", FALSE, NA, NA,
                    "Should caching of events or module be used?")
  ),
  inputObjects = bindrows(
    expectsInput(objectName = "cohortDataFiles", objectClass = "character",
                 desc = "A vector of file paths for cohortData output objects",
                 sourceURL = NA),
    expectsInput("ecoregionLayer", "SpatialPolygonsDataFrame",
                 desc = paste("A SpatialPolygonsDataFrame that characterizes the unique ecological regions used in the simulation",
                              "to parameterize the biomass, cover, and species establishment probability models. MUST be provided"),
                 sourceURL = NA),
    expectsInput(objectName = "pixelGroupMapFiles", objectClass = "character",
                 desc = "A vector of file paths for pixelGroupMap output objects",
                 sourceURL = NA),
    expectsInput(objectName = "rasterToMatch", objectClass = "character",
                 desc = paste("The rasterToMatch used in the simulation. MUST be provided"),
                 sourceURL = NA),
    expectsInput(objectName = "rstCurrentFiresFiles", objectClass = "character",
                 desc = "A vector of file paths for rstCurrentFires output objects",
                 sourceURL = NA),
    expectsInput(objectName = "severityDataFiles", objectClass = "character",
                 desc = "A vector of file paths for severityData output objects",
                 sourceURL = NA),
    expectsInput(objectName = "sppEquiv", objectClass = "character",
                 desc = paste("The sppEquiv used in the simulation. MUST be provided"),
                 sourceURL = NA),
    expectsInput(objectName = "vegTypeMapFiles", objectClass = "character",
                 desc = "A vector of file paths for vegTypeMap output objects",
                 sourceURL = NA)
  ),
  outputObjects = bindrows(
    createsOutput("allPixelBurnData", "data.table",
                  desc = paste("Pixelwise fire data per pixel: fire interval (no. years between each fire;",
                               "'fireInt'), fire frequency (mean `fireInt`; 'fireFreq'), severity (severity",
                               "class from *Biomass_regenerationPM*, assumed `5L` in stand-replacing ('noPM')",
                               "scenario), patch size (no. pixels with same severity class within a fire ID;",
                               "'patchSize') and killed total biomass ('severityB') and relative biomass",
                               "('severityPropB').")),
    createsOutput("allPixelCohortData", "data.table",
                  desc = paste("Pixelwise cohort data across the simulation landscape, plus total no fires",
                               "('noFires') and fire presence/absence ('firePresAbs').")),
    createsOutput("allPixelCohortDataMnt", "data.table",
                  desc = paste("As 'allPixelCohortData' but only on Montane region, with added vegetation",
                               "types from Cameron Naficy ('vegTypeCN')"))
  )
))

## event types
#   - type `init` is required for initialization

doEvent.LIM_resultsDataPrep = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      ### check for more detailed object dependencies:
      ### (use `checkObject` or similar)

      # do stuff for this event
      sim <- Init(sim)
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "loadVegData", eventPriority = 1)
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "loadFireData", eventPriority = 1.1)
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "calcFireMetrics", eventPriority = 2)
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "joinSimulationData", eventPriority = 3)
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "addVegTypesCN", eventPriority = 4)
    },
    loadVegData = {
      sim <- loadVegetationDataEvent(sim)
    },
    loadFireData = {
      sim <- loadFireDataEvent(sim)
    },
    calcFireMetrics = {
      sim <- calcFireAttributesEvent(sim)
    },
    joinSimulationData = {
      sim <- joinSimulationDataEvent(sim)
    },
    addVegTypesCN = {
      sim <- addVegTypesCNEvent(sim)
    },
    warning(paste("Undefined event type: \'", current(sim)[1, "eventType", with = FALSE],
                  "\' in module \'", current(sim)[1, "moduleName", with = FALSE], "\'", sep = ""))
  )
  return(invisible(sim))
}

### initialization
Init <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  mod$doAssertion <- getOption("LandR.assertions", TRUE)

  if (P(sim)$startYear > P(sim)$endYear) {
    stop("P(sim)$startYear can't be larger than P(sim)$endYear")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

loadVegetationDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "loadVegData")

  mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...

  grepPattrn <- paste0("/", P(sim)$scenarios)
  if (all(!is.na(P(sim)$reps))) {
    reps <- P(sim)$reps
    reps[reps < 10] <- paste0("0?", reps[reps < 10])
    grepPattrn <- paste0(grepPattrn, "_",
                         "rep", "(", paste0(reps, collapse = "|"), ")/")
  }
  names(grepPattrn) <- P(sim)$scenarios

  ## make list of stacked rasters -- filters between start and end years
  pixelGroupMapStkList <- Cache(Map,
                                x = grepPattrn,
                                f = loadStackFromRDS,
                                MoreArgs = list(files = sim$pixelGroupMapFiles,
                                                startYear = P(sim)$startYear,
                                                endYear = P(sim)$endYear),
                                userTags = c(cacheTags, "pixelGroupMapStkList"),
                                omitArgs = c("userTags"))

  vegTypeMapStkList <- Cache(Map,
                             x = grepPattrn,
                             f = loadStackFromRDS,
                             MoreArgs = list(files = sim$vegTypeMapFiles,
                                             startYear = P(sim)$startYear,
                                             endYear = P(sim)$endYear),
                             userTags = c(cacheTags, "vegTypeMapStkList"),
                             omitArgs = c("userTags"))


  ## pixelCohortData tables -- filters between start and end years and yearSubset
  cacheExtra <- sum(stack(sapply(pixelGroupMapStkList, function(ras) sum(ras))))
  cacheExtra <- sum(cacheExtra[], na.rm = TRUE)
  pixelCohortDataList <- Cache(Map,
                               x = grepPattrn,
                               f = loadCohortDataFromRDS,
                               pixelGroupMapStk = pixelGroupMapStkList,
                               MoreArgs = list(files = sim$cohortDataFiles,
                                               yearSubset = P(sim)$yearSubset,
                                               startYear = P(sim)$startYear,
                                               endYear = P(sim)$endYear),
                               .cacheExtra = cacheExtra,
                               userTags = c(cacheTags, "pixelCohortDataList"),
                               omitArgs = c("userTags", "pixelGroupMapStk"))

  ## vegTypeData tables -- filters to yearSubset
  cacheExtra2 <- sum(stack(sapply(vegTypeMapStkList, function(ras) sum(ras))))
  cacheExtra2 <- sum(cacheExtra[], na.rm = TRUE)
  vegTypeDataList <- Cache(Map,
                           f = vegTypeDataFromStks,
                           vegTypeMapStk = vegTypeMapStkList,
                           pixelGroupMapStk = pixelGroupMapStkList,
                           MoreArgs = list(yearSubset = P(sim)$yearSubset),
                           .cacheExtra = cacheExtra,
                           userTags = c(cacheTags, "vegTypeDataList"),
                           omitArgs = c("userTags", "vegTypeMapStk", "pixelGroupMapStk"))

  ## clean ws
  rm(vegTypeMapStkList, pixelGroupMapStkList)
  amc::.gc()

  ## join tables, add scenario col and rbind.
  ## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
  ## pixelBurntData can also have a different number of years if the saving frequency differs
  ## so join by keeping all pixels, calculate fire properties per pixel, then subset to
  ## cohort data years.
  pixelCohortDataList <- Cache(Map,
                               f = merge,
                               x = vegTypeDataList,
                               y = pixelCohortDataList,
                               MoreArgs = list(by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE),
                               .cacheExtra = c(cacheExtra, cacheExtra2),  ## use the same
                               userTags = c(cacheTags, "merge", "pixelCohortDataList"),
                               omitArgs = c("userTags", "x", "y"))
  amc::.gc()

  ## add scenario column when binding
  allPixelCohortData <- rbindlist(pixelCohortDataList, use.names = TRUE, fill = TRUE, idcol = "scenario")
  cols <- c("scenario",  "pixelIndex", "pixelGroup", "year", "rep", "vegType", "speciesCode",
            "ecoregionGroup", "age", "B", "mortality", "aNPPAct")   ## keep these columns only
  allPixelCohortData <- allPixelCohortData[, ..cols]
  if (exists("allPixelCohortData")) rm(pixelCohortDataList, vegTypeDataList)
  amc::.gc()

  ## export to sim
  sim$allPixelCohortData <- allPixelCohortData

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

loadFireDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "loadFireData")

  mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...

  grepPattrn <- paste0("/", P(sim)$scenarios)
  if (all(!is.na(P(sim)$reps))) {
    reps <- P(sim)$reps
    reps[reps < 10] <- paste0("0?", reps[reps < 10])
    grepPattrn <- paste0(grepPattrn, "_",
                         "rep", "(", paste0(reps, collapse = "|"), ")/")
  }
  names(grepPattrn) <- P(sim)$scenarios

  ## make list of stacked rasters -- filters between start and end years
  rstCurrentFiresStkList <- Cache(Map,
                                  x = grepPattrn,
                                  f = loadStackFromRDS,
                                  MoreArgs = list(files = sim$rstCurrentFiresFiles,
                                                  startYear = P(sim)$startYear,
                                                  endYear = P(sim)$endYear),
                                  userTags = c(cacheTags, "rstCurrentFiresStkList"),
                                  omitArgs = c("userTags"))

  ## severityData tables -- filters between start and end years
  severityDataList <- Cache(Map,
                            x = grepPattrn,
                            f = loadSeverityDataFromRDS,
                            MoreArgs = list(files = sim$severityDataFiles,
                                            startYear = P(sim)$startYear,
                                            endYear = P(sim)$endYear),
                            userTags = c(cacheTags, "severityDataList"),
                            omitArgs = c("userTags"))

  ## run function for each scenario then bind lists
  ## note that zeroes are added everywhere  where there is no severity data
  sevDataLs <- split(severityDataList$PM, by = c("year", "rep"))
  cacheExtra <- colSums(severityDataList$PM)
  severityRastersPM <- Cache(Map,
                             sevData = sevDataLs,
                             f = makeSevRasters,
                             MoreArgs = list(rasterToMatch = sim$rasterToMatch),
                             .cacheExtra = list(cacheExtra),
                             userTags = c(cacheTags, "severityRastersPM"),
                             omitArgs = c("userTags", "sevData"))

  ## add year and rep to names, to make the same as rstCurrentFiresStk
  names(severityRastersPM) <- paste0("year", sub("(.*)\\.(.*)", "\\1_rep\\2", names(severityRastersPM)))

  sevDataLs <- split(severityDataList$noPM, by = c("year", "rep"))
  cacheExtra <- colSums(severityDataList$noPM)
  severityRastersnoPM <- Cache(Map,
                               sevData = sevDataLs,
                               f = makeSevRasters,
                               MoreArgs = list(rasterToMatch = sim$rasterToMatch),
                               .cacheExtra = list(cacheExtra),
                               userTags = c(cacheTags, "severityRastersnoPM"),
                               omitArgs = c("userTags", "sevData"))
  names(severityRastersnoPM) <- paste0("year", sub("(.*)\\.(.*)", "\\1_rep\\2", names(severityRastersnoPM)))

  severityRasters <- list(PM = severityRastersPM, noPM = severityRastersnoPM)

  ## export to mod
  mod$severityRasters <- severityRasters
  mod$rstCurrentFiresStkList <- rstCurrentFiresStkList

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

calcFireAttributesEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  cacheTags <- c(currentModule(sim), "calcFireAttributes")

  ## FIRE ATTRIBUTES ---------------------------------------
  message(cyan("Calculating fire attributes..."))
  message(cyan("Patch size"))

  ## make rasters of patch size
  ## patchSizes are calculated only for fires that burned at least one forest pixel (otherwise there is no severity data)
  ## but consider the size of patches across forested and non-forested pixels (the severity of the later always being 0)
  ## then convert everything to a data.table for pixel level calculations
  ## by scenario and repetition
  cacheExtra <- sum(rast(lapply(mod$severityRasters$PM, function(ras) ras[[1]])), na.rm = TRUE)

  ## not all fire rasters have severity (in some years fires can remain outside forested pix), but
  ## all severity years should have a fire raster
  if (length(setdiff(names(mod$severityRasters$PM), names(mod$rstCurrentFiresStkList$PM)))) {
    stop("Not all fire severity rasters have an associated fire perimeter raster.")
  }

  # opts <- options(reproducible.useCache = FALSE)  ## test
  # on.exit(options(opts))

  ## for easier debugging
  # rasToDo <- names(mod$rstCurrentFiresStkList$PM)[1:3]   ## test
  rasToDo <- names(mod$rstCurrentFiresStkList$PM)
  missingRas <- setdiff(rasToDo, names(mod$severityRasters$PM))

  if (length(missingRas)) {
  ## add missing severity rasters from fires that did not burn forested pixels (severity is 0)
    tempRas <- mod$rstCurrentFiresStkList$PM[[missingRas]]
    if (inherits(mod$severityRasters$PM[[1]], "SpatRaster") &
        !inherits(tempRas, "SpatRaster")) {
      tempRas <- rast(tempRas)
    }
    tempRasLs <- lapply(tempRas, function(ras) {
      ras[!is.na(as.vector(ras[]))] <- 0
      rast(list(severityRas = ras, severityBRas = ras))
      })
    names(tempRasLs) <- names(tempRas)
    mod$severityRasters$PM <- c(mod$severityRasters$PM, tempRasLs)
  }

  tempList <- unstack(mod$rstCurrentFiresStkList$PM)   ## Map doesn't like to deal with different indexing of RasterStacks
  names(tempList) <- names(mod$rstCurrentFiresStkList$PM)

  patchSizeRasPM <- Cache(Map,
                          sevClassRasLs = mod$severityRasters$PM[rasToDo],
                          fireRas = tempList[rasToDo], ## subset and re-order to match
                          f = function(sevClassRasLs, fireRas) calcPatchSize(sevClassRasLs$severityRas, fireRas),
                          .cacheExtra = list(cacheExtra),
                          userTags = c(cacheTags, "patchSizeRasPM"),
                          omitArgs = c("userTags", "sevRasLs", "fireRasLs"))

  cacheExtra2 <- sum(rast(lapply(mod$severityRasters$noPM, function(ras) ras[[1]])))

  if (length(setdiff(names(mod$severityRasters$noPM), names(mod$rstCurrentFiresStkList$noPM)))) {
    stop("Not all fire severity rasters have an associated fire perimeter raster.")
  }

  ## for easier debugging
  # rasToDo <- names(mod$rstCurrentFiresStkList$noPM)[1:3]   ## test
  rasToDo <- names(mod$rstCurrentFiresStkList$noPM)
  missingRas <- setdiff(rasToDo, names(mod$severityRasters$noPM))

  if (length(missingRas)) {
    ## add missing severity rasters from fires that did not burn forested pixels (severity is 0)
    tempRas <- mod$rstCurrentFiresStkList$noPM[[missingRas]]
    if (inherits(mod$severityRasters$noPM[[1]], "SpatRaster") &
        !inherits(tempRas, "SpatRaster")) {
      tempRas <- rast(tempRas)
    }
    tempRasLs <- lapply(tempRas, function(ras) {
      ras[!is.na(as.vector(ras[]))] <- 0
      rast(list(severityRas = ras, severityBRas = ras))
    })
    names(tempRasLs) <- names(tempRas)
    mod$severityRasters$noPM <- c(mod$severityRasters$noPM, tempRasLs)
  }


  tempList <- unstack(mod$rstCurrentFiresStkList$noPM)   ## Map doesn't like to deal with different indexing of RasterStacks
  names(tempList) <- names(mod$rstCurrentFiresStkList$noPM)

  patchSizeRasnoPM <- Cache(Map,
                            sevClassRasLs = mod$severityRasters$noPM[rasToDo],
                            fireRas = tempList[rasToDo], ## subset and re-order to match
                            f = function(sevClassRasLs, fireRas) calcPatchSize(sevClassRasLs$severityRas, fireRas),
                            .cacheExtra = list(cacheExtra),
                            userTags = c(cacheTags, "patchSizeRasnoPM"),
                            omitArgs = c("userTags", "sevRasLs", "fireRasLs"))

  ## make a table of patch size -- use the same cacheExtra
  ## note that only pixels within fire perimeters are here (even if unforested and with 0 sev)
  patchSizeDataList <- Cache(Map,
                             fireAttrRasLs = list(PM = patchSizeRasPM, noPM = patchSizeRasnoPM),
                             f = fireAttrDTFromRasLs,
                             .cacheExtra = list(cacheExtra, cacheExtra2, "patchSize"),
                             userTags = c(cacheTags, "patchSizeDataList"),
                             omitArgs = c("userTags", "fireAttrRasLs"))
  allPatchSizeData <- rbindlist(patchSizeDataList, idcol = "scenario", fill = TRUE, use.names = TRUE)
  rm(patchSizeDataList, patchSizeRasnoPM, patchSizeRasPM, tempList); gc(reset = TRUE)


  ## make a table of fire severity -- use the same cacheExtra
  message(cyan("Severity"))
  ## note that all pixels are here, burnt and unburnt
  ## for easier debugging
  # rasToDo <- list(PM = names(mod$severityRasters$PM)[1:3],
  #                 noPM = names(mod$severityRasters$noPM)[1:3]) ## test
  rasToDo <- list(PM = names(mod$severityRasters$PM),
                  noPM = names(mod$severityRasters$noPM))
  rasToDo <- rasToDo[names(mod$severityRasters)] ## ensure order is correct
  severityDataList <- Cache(Map,
                            fireAttrRasLs = mod$severityRasters,
                            i = rasToDo,
                            f = fireAttrDTFromRasLs,
                            .cacheExtra = list(cacheExtra, cacheExtra2, "severity"),
                            userTags = c(cacheTags, "severityDataList"),
                            omitArgs = c("userTags", "fireAttrRasLs"))
  ## add scenario column when binding
  allSeverityData <- rbindlist(severityDataList, idcol = "scenario", fill = TRUE, use.names = TRUE)
  setnames(allSeverityData, c("severityRas", "severityBRas"), c("severity", "severityB"))

  rm(severityDataList); gc(reset = TRUE)

  ## make table with fire occurrences and IDs per year/rep -- use same cacheExtra
  message(cyan("Fire occurrences"))

  ## for easier debugging
  # rasToDo <- list(PM = names(mod$rstCurrentFiresStkList$PM)[1:3],
  #                 noPM = names(mod$rstCurrentFiresStkList$noPM)[1:3]) ## testing
  rasToDo <- list(PM = names(mod$rstCurrentFiresStkList$PM),
                  noPM = names(mod$rstCurrentFiresStkList$noPM))
  rasToDo <- rasToDo[names(mod$rstCurrentFiresStkList)] ## to ensure list order is the same
  pixelBurnDataList <- Cache(Map,
                             rstCurrentFiresStk = mod$rstCurrentFiresStkList,
                             i = rasToDo,
                             f = pixelBurnDataFromStks,
                             .cacheExtra = list(cacheExtra, cacheExtra2),
                             userTags = c(cacheTags, "pixelBurnDataList"),
                             omitArgs = c("userTags", "rstCurrentFiresStk"))
  ## add scenario column when binding
  allPixelBurnData <- rbindlist(pixelBurnDataList, idcol = "scenario", use.names = TRUE)
  rm(pixelBurnDataList); gc(reset = TRUE)

  if (mod$doAssertion)  {
    test <- allPixelBurnData[!allPatchSizeData, on = .(scenario, rep, year, pixelIndex)]  ## should be empty of burnt pixels.
    if (nrow(test)) {
      stop("Some burnt pixels (per scenario/rep/year) have no patch sizes associated to them")
    }

    test <- allPatchSizeData[!allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]  ## should be empty of burnt pixels.
    if (nrow(test)) {
      stop("Some pixels with patch sizes (per scenario/rep/year) have no fire ID")
    }
    rm(test); gc(reset = TRUE)
  }

  ## join other tables
  allPixelBurnData <- allPatchSizeData[allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]
  rm(allPatchSizeData); gc(reset = TRUE)   ## clean progressively to release memory

  ## keep only the pixels that burned here -- saves memory if we don't have non-burns per year
  ## which are NAs in mean severity/patch size calculations.
  allPixelBurnData <- allPixelBurnData[allSeverityData, nomatch = 0, on = .(scenario, rep, year, pixelIndex)]
  rm(allSeverityData); gc(reset = TRUE)

  ## no. fires per pixel
  ## how many times did each pixel burn? total no. fires per pixel/scenario/rep
  allPixelBurnData[, noFires := as.integer(sum(burnt, na.rm = TRUE)), by = .(scenario, rep, pixelIndex)]

  ## calculate fire size in pixels per fireID/scenario/rep
  allPixelBurnData[!is.na(fireID), patchSizePix := as.integer(length(unique(pixelIndex))),
                    by = .(scenario, rep, year, fireID)]
  ## rename log_area to patchSizeLogHa
  setnames(allPixelBurnData, "log_area", "patchSizeLogHa")

  ## Fire intervals -------
  ## calculate fire frequency as mean fire intervals per scenario, rep, pixel, considering start and end years
  ## (mean FRI in Steel et al 2021, see this paper for limitations and details)
  message(cyan("Fire intervals"))

  ## calculate intervals to a new table (there will be always one interval more than fires)
  ## first calculate intervals in pixels with a fire history
  fireIntervals <-  allPixelBurnData[noFires > 0,
                                     list(fireInt = as.integer(c(P(sim)$startYear, year, P(sim)$endYear) - lag(c(P(sim)$startYear, year, P(sim)$endYear), n = 1)),
                                          noFires = unique(noFires)),
                                     by = .(scenario, rep, pixelIndex)]
  fireIntervals <- fireIntervals[!is.na(fireInt)]   ## NAs come from startYear, and can be removed
  ## check
  if (mod$doAssertion) {
    if (any(fireIntervals[, length(fireInt) != unique(noFires) + 1, by = .(scenario, rep, pixelIndex)]$V1))
      stop("Fire interval calculations are wrong")
  }

  ## fire frequency == average fire return interval as in Steel et al 2021
  fireIntervals[, fireFreq := mean(fireInt, na.rm = TRUE), by = .(scenario, rep, pixelIndex)]
  fireIntervals <- unique(fireIntervals[, .(scenario, rep, pixelIndex, fireFreq)])

  ## join to datatable
  setkey(fireIntervals, pixelIndex, scenario, rep)
  setkey(allPixelBurnData, pixelIndex, scenario, rep)

  allPixelBurnData <- fireIntervals[allPixelBurnData]

  rm(fireIntervals)
  amc::.gc()

  ## PIXELS WITH NO FIRE HISTORY ------------
  ## make a table of pixels that have never burned, per scenario
  noFireHistoryDataLsPM <- Cache(Map,
                                 r = P(sim)$reps,
                                 MoreArgs = list(
                                   rstCurrentFiresStk = mod$rstCurrentFiresStkList$PM,
                                   rasterToMatch = sim$rasterToMatch,
                                   doAssertion = mod$doAssertion),
                                 f = makeNoFireHistoryData,
                                 .cacheExtra = list(cacheExtra, cacheExtra2),
                                 userTags = c(cacheTags, "noFireHistoryDataPM"),
                                 omitArgs = c("userTags", "rstCurrentFiresStk"))
  noFireHistoryDataPM <- rbindlist(noFireHistoryDataLsPM, use.names = TRUE)
  rm(noFireHistoryDataLsPM); gc(reset = TRUE)

  noFireHistoryDataLsnoPM <- Cache(Map,
                                 r = P(sim)$reps,
                                 MoreArgs = list(
                                   rstCurrentFiresStk = mod$rstCurrentFiresStkList$noPM,
                                   rasterToMatch = sim$rasterToMatch,
                                   doAssertion = mod$doAssertion),
                                 f = makeNoFireHistoryData,
                                 .cacheExtra = list(cacheExtra, cacheExtra2),
                                 userTags = c(cacheTags, "noFireHistoryDatanoPM"),
                                 omitArgs = c("userTags", "rstCurrentFiresStk"))
  noFireHistoryDatanoPM <- rbindlist(noFireHistoryDataLsnoPM, use.names = TRUE)
  rm(noFireHistoryDataLsnoPM); gc(reset = TRUE)

  ## TODO: LEFT OFF HERE.
  noFireHistoryData <- rbindlist(list(PM = noFireHistoryDataPM, noPM = noFireHistoryDatanoPM),
                                 use.names = TRUE, idcol = "scenario")
  rm(noFireHistoryDataPM, noFireHistoryDatanoPM); gc(reset = TRUE)

  ## add other fire attributes for pixels with no fire history
  ## make column of noFires
  noFireHistoryData[, `:=`(noFires = 0L, fireFreq = P(sim)$endYear - P(sim)$startYear,
                           severity = 0L, severityB = 0, burnt = 0L)]
  noFireHistoryData[, ID := NULL]

  ## check
  if (mod$doAssertion) {
    test <- allPixelBurnData[noFireHistoryData, nomatch = 0, on = .(scenario, rep, pixelIndex)]
    if (nrow(test)) {
      stop("There should be no common scenario/rep/pixelIndex combinations between fire history data and no-fire-history data")
    }
    rm(test); gc(reset = TRUE)
  }

  ## bind pixels that never burned
  allPixelBurnData <- rbindlist(list(allPixelBurnData, noFireHistoryData), fill = TRUE, use.names = TRUE)

  # OLD CODE
  # ## calculate patch size, as the number of in pixels per severity (class)/fireID/scenario/rep
  # ## note that for noPM we assume severity class (i.e. 'severity' column) to be the maximum = 5
  # ## only in pixels with a pixelGroup (others are non-forest and had no veg dynamics)
  # ## also note that we are ignoring if patches are contiguous or not, and simply counting the number of pixels
  # ## with a given severity per fireID
  # message(blue("Assuming a severity class 5 for any scenario with 'noPM'"))
  # allPixelBurnData[grepl("noPM", scenario) & !is.na(pixelGroup), severity := 5]
  # allPixelBurnData[, severity := as.integer(severity)]
  # allPixelBurnData[!is.na(severity), patchSize := as.integer(length(unique(pixelIndex))),
  #                  by = .(scenario, rep, year, severity, fireID)]

  ## fire frequency
  ## calculate fire frequency as the mean fire-intervals per pixel (see Steel et al 2021 for limitations and details)
  # setkey(allPixelBurnData, pixelIndex, scenario, rep, year)
  # allPixelBurnData[, fireInt := as.integer(year - lag(year, n = 1)),
  #                  by = .(scenario, rep, pixelIndex)]
  # allPixelBurnData[is.na(fireInt), fireInt := as.integer(year - P(sim)$startYear)] ## NAs mean only one fire, return interval is the difference from start year
  # allPixelBurnData[, fireFreq := mean(fireInt), by = .(scenario, rep, pixelIndex)]
  #
  # allPixelBurnData[, burnt := NULL] ## no longer necessary

  ## checks
  if (mod$doAssertion)  {
    test1 <- sapply(split(allPixelBurnData, by = c("scenario", "rep", "year")), FUN = function(x){
      any(duplicated(x[, pixelIndex]))
    })
    if (any(test1))
      stop("Each pixel should only have one record of no. fires per scenario")

    ## OLD CODE
    # test2 <- setdiff(which(is.na(allPixelBurnData$pixelGroup)),
    #                  which(is.na(allPixelBurnData$severity)))
    # test3 <- setdiff(which(is.na(allPixelBurnData$pixelGroup)),
    #                  which(is.na(allPixelBurnData$severityB)))
    # if (length(test2) | length(test3)) {
    #   stop("NAs differ between pixelGroup and severity/severityB")
    # }

    test4 <- any(is.na(allPixelBurnData[!is.na(severity), patchSizeLogHa])) |
      any(is.na(allPixelBurnData[!is.na(severity), patchSizePix]))
    if (test4) {
      stop("There are NA's in patch sizes where severity is non-NA")
    }

    test5 <- any(is.na(allPixelBurnData$fireFreq))
    if (test5) {
      stop("NA fire intervals")
    }
    suppressWarnings(rm(test1, test2, test3, test4, test5))
  }

  sim$allPixelBurnData <- allPixelBurnData

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

joinSimulationDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  gc(reset = TRUE)  ## try to release memory consumed by DT threads
  mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...

  ## ECOLOGICAL ZONATION -----------------------------
  message(cyan("Adding ecological zones"))
  ecoregionLayerRas <- rasterize(sim$ecoregionLayer, sim$rasterToMatch, field = "ecozoneCode")
  ecoregionLayerDT <- data.table(ecozoneCode = getValues(ecoregionLayerRas),
                                 pixelIndex = seq_len(ncell(ecoregionLayerRas)))
  ecoregionLayerDT <- ecoregionLayerDT[!is.na(ecozoneCode)]

  ecoregionLayerLabels <- data.table(ecozoneCode = sim$ecoregionLayer$ecozoneCode,
                                     ecozoneName = paste(sim$ecoregionLayer$NRNAME,
                                                         sim$ecoregionLayer$NSRNAME, sep = " - ")) %>%
    unique(.)

  ecoregionLayerDT <- ecoregionLayerLabels[ecoregionLayerDT, on = .(ecozoneCode)]

  allPixelCohortData <- ecoregionLayerDT[sim$allPixelCohortData, on = .(pixelIndex)]
  allPixelCohortData[, `:=`(scenario = as.factor(scenario),
                            rep = as.integer(rep),
                            year = as.integer(year),
                            pixelGroup = as.integer(pixelGroup),
                            pixelIndex = as.integer(pixelIndex),
                            ecoregionGroup = as.factor(ecoregionGroup),
                            ecozoneCode = as.integer(ecozoneCode),
                            ecozoneName = as.factor(ecozoneName),
                            vegType = as.integer(vegType),
                            speciesCode = as.factor(speciesCode),
                            age = as.integer(age),
                            B = as.integer(B),
                            mortality = as.integer(mortality),
                            aNPPAct = as.integer(aNPPAct))]
  amc::.gc()

  ## add noFires to cohortData - no year info, because its the total across the simulation
  cols <- c("scenario", "rep", "pixelIndex", "noFires")
  allPixelCohortData <- unique(sim$allPixelBurnData[, ..cols])[allPixelCohortData,
                                                               on = .(scenario, rep, pixelIndex)]

  ## checks
  if (mod$doAssertion) {
    cols <- c("scenario", "rep", "pixelIndex")
    test1 <- sim$allPixelBurnData[allPixelCohortData[is.na(noFires), ..cols], on = cols, nomatch = 0]
    test1 <- dim(test1)[1]
    if (test1) {
      stop("There shouldn't be any NAs in noFires except in scenario/rep/pixelIndex\n",
           "combos that did not have any fires during the simulation")
    }
    rm(test1)
  }

  allPixelCohortData[is.na(noFires), noFires := 0]
  allPixelCohortData <- allPixelCohortData[!is.na(pixelGroup),]
  amc::.gc()

  ## ADD MISSING SPECIES IN YEAR/SCENARIO/PIXEL COMBINATION
  ## cohortData doens't track absent cohorts, so they need to ba added back
  ## for now pixels from the 0s pixelGroup  have one entry with NAs for speciesCode
  ## they will be ignored for now and removed later, after adding one species entry for each of these pixels.
  ## for reporting consistency add to show losses in B
  message(cyan("Adding absent species in all scenario/year/pixel combinations"))
  combinations <- unique(allPixelCohortData[, .(scenario, rep, year, pixelIndex, pixelGroup)])
  spp <- as.character(na.omit(unique(allPixelCohortData$speciesCode)))
  combinations <- lapply(spp, FUN = function(x) {
    data.table(combinations,
               speciesCode = x)
  }) %>%
    rbindlist(., use.names = TRUE)

  ## join while keeping all combos, NA species will now disappear.
  allPixelCohortData <- allPixelCohortData[combinations,
                                           on = .(scenario, rep, year, pixelIndex,
                                                  pixelGroup, speciesCode)]
  rm(spp, combinations)
  amc::.gc()

  ## checks
  if (mod$doAssertion) {
    test <- length(unique(allPixelCohortData[, length(unique(pixelIndex)), by = .(scenario, rep, year)]$V1)) == 1
    test2 <- length(unique(allPixelCohortData[, length(unique(speciesCode)), by = .(scenario, rep, year, pixelIndex)]$V1)) == 1
    test3 <- any(is.na(allPixelCohortData$speciesCode))

    if (isFALSE(test))
      stop("No. pixels should be the same across years, for a given scenario/rep")
    if (isFALSE(test2))
      stop("No. species per pixel should be the same across pixels, for a given scenario/rep/year")
    if (test3)
      stop("There are NA speciesCodes")
    rm(test, test2, test3)
    amc::.gc()
  }

  ## replace NAs of cohortData by 0s
  cols <- c("age", "B", "mortality", "aNPPAct", "vegType", "noFires")
  allPixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
  amc::.gc()

  ## add ecoregion group/ecozone code/name where they're missing
  ## add vegType where it's missing, but it's a pixel with some veg
  allPixelCohortData[, `:=`(vegType = max(vegType)),
                     by = .(scenario, rep, year, pixelGroup)]

  allPixelCohortData[, `:=`(ecoregionGroup = unique(na.omit(ecoregionGroup)),
                            ecozoneCode = unique(na.omit(ecozoneCode)),
                            ecozoneName = unique(na.omit(ecozoneName))),
                     by = .(pixelIndex)]
  amc::.gc()

  ## add noFires where it's missing
  allPixelCohortData[, noFires := max(noFires),
                     by = .(scenario, rep, pixelIndex)]
  amc::.gc()

  ## add presence/absence of fire across simulation per pixel/scenario
  allPixelCohortData[, firePresAbs := as.integer(any(noFires > 0)),
                     by = .(scenario, rep, pixelIndex)]
  amc::.gc()

  ## export to sim
  sim$allPixelCohortData <- allPixelCohortData

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

addVegTypesCNEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  gc()  ## try to release memory consumed by DT threads
  mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...

  ## USING CAMERON'S CLASSIFICATION/SUMMARY ---------------------
  ## Cameron uses relative basal area to classify stand structure, we can use relative Biomass.
  ## we subset to the montane ecological zone, from where Cameron's data comes from
  message(cyan("Classifying vegetation types according to Cameron's classification"))

  allPixelCohortDataMnt <- sim$allPixelCohortData[grep("Montane", ecozoneName)]
  allPixelCohortDataMnt[, sumB := sum(B), by = .(scenario, year, rep, pixelGroup)]
  allPixelCohortDataMnt[, relB := sum(B)/sumB, by = .(scenario, rep, year, pixelGroup, speciesCode)]
  allPixelCohortDataMnt[is.na(relB) & sumB == 0, relB := 0]
  amc::.gc()

  if (any(is.na(allPixelCohortDataMnt$relB)))
    stop("Missing values in relative biomass")

  ## subset to a smaller DT and join Cameron's species names
  vegTypesCN <- unique(allPixelCohortDataMnt[B > 0, .(scenario, rep, year, pixelGroup, speciesCode, relB)])
  vegTypesCN <- unique(na.omit(sim$sppEquiv[, .(Cameron, LIM)]))[vegTypesCN, on = "LIM==speciesCode",
                                                                 allow.cartesian = TRUE]
  setnames(vegTypesCN, "LIM", "speciesCode")

  amc::.gc()
  if (P(sim)$parallel) {
    if (.Platform$OS.type == "windows") {
      plan("multisession", workers = P(sim)$ncores)
    } else {
      plan("multicore", workers = P(sim)$ncores)
    }
    vegTypesCN <- future_lapply(split(vegTypesCN, by = c("scenario", "rep")),
                                FUN = parallelFUN, cPath = cachePath(sim))
    future:::ClusterRegistry("stop")
  } else {
    vegTypesCN <- lapply(split(vegTypesCN, by = c("scenario", "rep")),
                         FUN = parallelFUN, cPath = cachePath(sim))
  }

  amc::.gc()
  vegTypesCN <- rbindlist(vegTypesCN, use.names = TRUE)

  ## test:
  # vegTypesCN <- lapply(unique(showCache(cachePath(sim),
  #                                       userTags = c("convertToCNVegType", "reportResults"),
  #                                       after = "2021-06-16")$cacheId),
  #                      FUN = function(x) {loadFromCache(cachePath(sim), cacheId = x)}) %>%
  #   rbindlist(.)

  ## add Cameron's veg types and get rid of useless columns
  cols <- c("scenario", "rep", "year", "pixelGroup", "vegTypeCN")
  cols2 <- c("scenario", "rep", "year", "pixelGroup")
  allPixelCohortDataMnt <- tryCatch(unique(vegTypesCN[, ..cols])[allPixelCohortDataMnt,
                                                                 on = cols2],
                                    error = allPixelCohortDataMnt)
  if (!"vegTypeCN" %in% names(allPixelCohortDataMnt)) {
    stop("Joining Cameron's veg types didn't work")
  }

  if (any(is.na(allPixelCohortDataMnt$vegTypeCN) & allPixelCohortDataMnt$B > 0)) {
    stop("Some pixels with biomass were not assigned a vegTypeCN")
  }
  rm(vegTypesCN)

  allPixelCohortDataMnt[is.na(vegTypeCN), vegTypeCN := "No veg."]
  allPixelCohortDataMnt[, `:=`(sumB = NULL,
                               relB = NULL,
                               vegType = NULL)]
  ## make "No veg." the last factor
  levs <- c(sort(grep("No veg.", unique(allPixelCohortDataMnt$vegTypeCN), value = TRUE, invert = TRUE)),
            "No veg.")
  allPixelCohortDataMnt[, vegTypeCN := factor(vegTypeCN, levels = levs)]
  amc::.gc()

  ## export to sim
  allPixelCohortDataMnt[, `:=`(scenario = as.factor(scenario),
                               rep = as.integer(rep),
                               year = as.integer(year),
                               pixelGroup = as.integer(pixelGroup),
                               pixelIndex = as.integer(pixelIndex),
                               ecoregionGroup = as.factor(ecoregionGroup),
                               ecozoneCode = as.integer(ecozoneCode),
                               ecozoneName = as.factor(ecozoneName),
                               speciesCode = as.factor(speciesCode),
                               vegTypeCN = as.factor(vegTypeCN),
                               noFires = as.integer(noFires),
                               firePresAbs = as.integer(firePresAbs),
                               age = as.integer(age),
                               B = as.integer(B),
                               mortality = as.integer(mortality),
                               aNPPAct = as.integer(aNPPAct))]
  sim$allPixelCohortDataMnt <- allPixelCohortDataMnt

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

.inputObjects <- function(sim) {
  # Any code written here will be run during the simInit for the purpose of creating
  # any objects required by this module and identified in the inputObjects element of defineModule.
  # This is useful if there is something required before simulation to produce the module
  # object dependencies, including such things as downloading default datasets, e.g.,
  # downloadData("LCC2005", modulePath(sim)).
  # Nothing should be created here that does not create a named object in inputObjects.
  # Any other initiation procedures should be put in "init" eventType of the doEvent function.
  # Note: the module developer can check if an object is 'suppliedElsewhere' to
  # selectively skip unnecessary steps because the user has provided those inputObjects in the
  # simInit call, or another module will supply or has supplied it. e.g.,
  # if (!suppliedElsewhere('defaultColor', sim)) {
  #   sim$map <- Cache(prepInputs, extractURL('map')) # download, extract, load file from url in sourceURL
  # }

  #cacheTags <- c(currentModule(sim), "function:.inputObjects") ## uncomment this if Cache is being used
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # ! ----- EDIT BELOW ----- ! #
  grepPattrn <- paste0(P(sim)$scenarios, collapse = "|")
  if (all(!is.na(P(sim)$reps))) {
    reps <- P(sim)$reps
    reps[reps < 10] <- paste0("0?", reps[reps < 10])
    grepPattrn <- paste0("/(", grepPattrn, ")", "_",
                         "rep", "(", paste0(reps, collapse = "|"), ")$")
  } else {
    if (length(P(sim)$reps) > 1) {
      stop("P(sim)$reps has NAs and a length > 1. Please supply:\n",
           "P(sim)$reps <- NA OR a non-NA integer vector")
    }
  }
  outputDirs <- list.dirs(outputPath(sim), full.names = TRUE, recursive = TRUE) %>%
    grep(grepPattrn, ., value = TRUE)

  ## GET FILE NAMES IF NOT SUPPLIED
  if (!suppliedElsewhere("cohortDataFiles", sim)) {
    sim$cohortDataFiles <- lapply(outputDirs, FUN = function(x) list.files(x, pattern = "cohortData", full.names = TRUE)) %>%
      do.call(c,.) %>%
      unique(.)
  }

  if (!suppliedElsewhere("severityDataFiles", sim)) {
    sim$severityDataFiles <- lapply(outputDirs, FUN = function(x) list.files(x, pattern = "severityData", full.names = TRUE)) %>%
      do.call(c,.) %>%
      unique(.)
  }

  if (!suppliedElsewhere("rstCurrentFiresFiles", sim)) {
    sim$rstCurrentFiresFiles <- lapply(outputDirs, FUN = function(x) list.files(x, pattern = "rstCurrentFires", full.names = TRUE)) %>%
      do.call(c,.) %>%
      unique(.)
  }

  if (!suppliedElsewhere("pixelGroupMapFiles", sim)) {
    sim$pixelGroupMapFiles <- lapply(outputDirs, FUN = function(x) list.files(x, pattern = "pixelGroupMap", full.names = TRUE)) %>%
      do.call(c,.) %>%
      unique(.)
  }

  if (!suppliedElsewhere("vegTypeMapFiles", sim)) {
    sim$vegTypeMapFiles <- lapply(outputDirs, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE)) %>%
      do.call(c,.) %>%
      unique(.)
  }

  if (!suppliedElsewhere("ecoregionLayer", sim)) {
    stop("Please provide 'ecoregionLayer' used in the simulation")
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    stop("Please provide 'sppEquiv' used in the simulation")
  }

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    stop("Please provide 'rasterToMatch' used in the simulation")
  }

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
