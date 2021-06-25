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
  version = list(SpaDES.core = "1.0.8.9000",
                 LIM_resultsDataPrep = "0.0.0.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = deparse(list("README.md", "LIM_resultsDataPrep.Rmd")), ## same file
  reqdPkgs = list("data.table", "raster", "LandR",
                  "future", "future.apply",
                  "CeresBarros/ToolsCB", "crayon"),
  parameters = rbind(
    defineParameter("endYear", "integer", 2111L, 1, NA,
                    "The last year of the simulation"),
    defineParameter("ncores", "integer", 8, 1, NA,
                    "Number of cores to use if P(sim)$parallel is TRUE"),
    defineParameter("parallel", "logical", TRUE, 1, NA,
                    paste("Should data processing be parallelized? Currently only used in assigning",
                          "functional vegetation types following Cameron Naficy's classification")),
    defineParameter("reps", "integer", 1L:10L, NA, NA,
                    "The simulation repetitions to compile. If no repetitions were performed set to NA"),
    defineParameter("scenarios", "character", c("PM", "noPM"), NA, NA,
                    paste("The simulation scenarios to compile - must correspond with the names used in the",
                          "outputs folder tree")),
    defineParameter("startYear", "integer", 2011L, NA, NA,
                    "The first year of the simulation"),
    defineParameter("yearSubset", "integer", unique(c(seq(2011L, 2111L, 5), 2111L)), NA, NA,
                    paste("The simulation years to compile - only vegetation dynamics will be subset.",
                          "outputs folder tree. If using all years set to NULL")),
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
    createsOutput("allPixelBurnData", "data.table", "Pixelwise fire data across the simulation landscape."),
    createsOutput("allPixelCohortData", "data.table", "Pixelwise cohort and fire data across the simulation landscape."),
    createsOutput("allPixelCohortDataMnt", "data.table", "Pixelwise cohort and fire data for montane belt.")
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
                           eventType = "loadSimulationData")
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "joinSimulationData")
      sim <- scheduleEvent(sim, eventTime = start(sim), moduleName = currentModule(sim),
                           eventType = "addVegTypesCN")
    },
    loadSimulationData = {
      sim <- loadSimulationDataEvent(sim)
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

loadSimulationDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  mod$doAssertion <- getOption("LandR.assertions", TRUE)  ## this is not being cached...

  grepPattrn <- paste0("/", P(sim)$scenarios)
  if (all(!is.na(P(sim)$reps))) {
    reps <- P(sim)$reps
    reps[reps < 10] <- paste0("0", reps[reps < 10])
    grepPattrn <- paste0(grepPattrn, "_",
                         "rep", "(", paste0(reps, collapse = "|"), ")/")
  }
  names(grepPattrn) <- P(sim)$scenarios

  rstCurrentFiresStkList <- sapply(grepPattrn,
                                   FUN = loadStackFromRDS,
                                   files = sim$rstCurrentFiresFiles,
                                   simplify = FALSE, USE.NAMES = TRUE)

  pixelGroupMapStkList <- sapply(grepPattrn,
                                 FUN = loadStackFromRDS,
                                 files = sim$pixelGroupMapFiles,
                                 simplify = FALSE, USE.NAMES = TRUE)

  vegTypeMapStkList <- sapply(grepPattrn,
                              FUN = loadStackFromRDS,
                              files = sim$vegTypeMapFiles,
                              simplify = FALSE, USE.NAMES = TRUE)


  ## pixelCohortData tables
  pixelCohortDataList <- mapply(FUN = loadCohortDataFromRDS,
                                x = grepPattrn,
                                pixelGroupMapStk = pixelGroupMapStkList,
                                MoreArgs = list(files = sim$cohortDataFiles,
                                                yearSubset = P(sim)$yearSubset),
                                SIMPLIFY = FALSE, USE.NAMES = TRUE)

  ## vegTypeData tables
  vegTypeDataList <- mapply(FUN = vegTypeDataFromStks,
                            vegTypeMapStk = vegTypeMapStkList,
                            pixelGroupMapStk = pixelGroupMapStkList,
                            MoreArgs = list(yearSubset = P(sim)$yearSubset),
                            SIMPLIFY = FALSE, USE.NAMES = TRUE)

  ## pixelBurnData tables - all rasters
  pixelBurnDataList <- sapply(rstCurrentFiresStkList,
                              FUN = pixelBurnDataFromStks,
                              simplify = FALSE, USE.NAMES = TRUE)
  ## add scenario column when binding
  allPixelBurnData <- rbindlist(pixelBurnDataList, idcol = "scenario", use.names = TRUE)
  allPixelBurnData <- allPixelBurnData[fireID != "NA"]

  ## severityData tables
  severityDataList <- sapply(grepPattrn,
                             FUN = loadSeverityDataFromRDS,
                             files = sim$severityDataFiles,
                             simplify = FALSE, USE.NAMES = TRUE)
  ## add scenario column when binding
  allSeverityData <- rbindlist(severityDataList, idcol = "scenario", fill = TRUE, use.names = TRUE)

  ## join fire data
  amc::.gc()
  allPixelBurnData <- allSeverityData[allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]

  ## clean ws - keep for labels
  rm(rstCurrentFiresStkList, vegTypeMapStkList, pixelGroupMapStkList,
     pixelBurnDataList, severityDataList, allSeverityData)
  amc::.gc()

  ## join tables, add scenario col and rbind.
  ## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
  ## pixelBurntData can also have a different number of years if the saving frequency differs
  ## so join by keeping all pixels, calculate fire properties per pixel, then subset to
  ## cohort data years.
  pixelCohortDataList <- mapply(FUN = merge,
                                x = vegTypeDataList,
                                y = pixelCohortDataList,
                                MoreArgs = list(by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE),
                                SIMPLIFY = FALSE, USE.NAMES = TRUE)
  amc::.gc()

  ## add scenario column when binding
  allPixelCohortData <- rbindlist(pixelCohortDataList, use.names = TRUE, idcol = "scenario")
  if (exists("allPixelCohortData")) rm(pixelCohortDataList, vegTypeDataList)
  amc::.gc()

  ## checks
  if (mod$doAssertion)  {
    ## note that some years may not have all the reps in allPixelBurnData if there where no fires
    ## because we excluded NAs (i.e. pixels without fires)
    repsFire <- unique(allPixelBurnData[, length(unique(rep)), by = .(scenario)]$V1)
    repsCohortData <- unique(allPixelCohortData[, length(unique(rep)), by = .(scenario, year)]$V1)
    test1 <- all(length(repsFire) == 1, length(repsCohortData) == 1)
    test2 <- all(length(repsFire) == length(repsCohortData),
                 identical(sort(repsFire), sort(repsCohortData)))

    if (!isTRUE(test1)) {
      stop("Fire and/or cohort data do not have the same number of reps across scenarios")
    }
    if (!isTRUE(test2)) {
      stop("Fire and cohort data differ in number of reps (per scenario/year)")
    }
  }

  ## export to sim
  sim$allPixelBurnData <- allPixelBurnData
  sim$allPixelCohortData <- allPixelCohortData

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

joinSimulationDataEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
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
  amc::.gc()

  ## FIRE ATTRIBUTES ---------------------------------------
  message(cyan("Calculating and adding fire attributes"))
  ## no. fires per pixel
  ## how many times did each pixel burn? total no. fires per pixel/scenario/rep
  allPixelBurnData <- copy(sim$allPixelBurnData)
  allPixelBurnData[, noFires := sum(burnt), by = .(scenario, rep, pixelIndex)]

  ## calculate fire size in pixels per fireID/scenario/rep
  ## this accounts for both forest and non-forest pixels
  allPixelBurnData[, fireSize := length(unique(pixelIndex)), by = .(scenario, rep, year, fireID)]

  ## calculate patch size, as the number of in pixels per severity (class)/fireID/scenario/rep
  ## note that for noPM we assume severity class (i.e. 'severity' column) to be the maximum = 5
  ## only in pixels with a pixelGroup (others are non-forest and had no veg dynamics)
  ## also note that we are ignoring if patches are contiguous or not, and simply counting the number of pixels
  ## with a given severity per fireID
  message(blue("Assuming a severity class 5 for any scenario with 'noPM'"))
  allPixelBurnData[grepl("noPM", scenario) & !is.na(pixelGroup), severity := 5]
  allPixelBurnData[!is.na(severity), patchSize := length(unique(pixelIndex)),
                   by = .(scenario, rep, year, severity, fireID)]

  ## fire frequency
  ## calculate fire frequency as the mean fire-intervals per pixel (see Steel et al 2021 for limitations and details)
  fireFreqDT <- allPixelBurnData[, list(year = c(year, P(sim)$endYear),     ## year one is dropped here
                                        fireInt = diff(c(P(sim)$startYear, year, P(sim)$endYear))),   ## interval calculated between fire years, and start and end years
                                 by = .(scenario, rep, pixelIndex)]
  ## because we forced a start and end year, intervals of 0 for the 100th year mean that there was a fire at year P(sim)$endYear
  ## this doesn't apply in the same way to fires at P(sim)$startYear, these should have a return interval of 0 (because we added the first year)
  ## if only one fire occurred and it was at P(sim)$endYear, then the correct interval is 99 (P(sim)$endYear-P(sim)$startYear)
  fireFreqDT <- fireFreqDT[!(fireInt == 0 & year == P(sim)$endYear)]
  fireFreqDT[, fireFreq := mean(fireInt), by = .(scenario, rep, pixelIndex)]

  ## join DTs
  allPixelBurnData <- fireFreqDT[allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]
  allPixelBurnData[, burnt := NULL] ## no longer necessary

  ## checks
  if (mod$doAssertion)  {
    test1 <- sapply(split(allPixelBurnData, by = c("scenario", "rep", "year")), FUN = function(x){
      any(duplicated(x[, pixelIndex]))
    })
    test2 <- setdiff(which(is.na(allPixelBurnData$pixelGroup)),
                     which(is.na(allPixelBurnData$severity)))
    test3 <- setdiff(which(is.na(allPixelBurnData$pixelGroup)),
                     which(is.na(allPixelBurnData$severityB)))
    test4 <- any(is.na(allPixelBurnData[!is.na(severity), patchSize]))
    test5 <- any(is.na(allPixelBurnData$fireFreq))

    if (any(test1))
      stop("Each pixel should only have one record of no. fires per scenario")

    if (length(test2) | length(test3)) {
      stop("NAs differ between pixelGroup and severity/severityB")
    }

    if (test4) {
      stop("There are NA's in patch sizes where severity is non-NA")
    }

    if (test5) {
      stop("NA fire intervals")
    }
    rm(test1, test2, test3, test4, test5)
  }
  rm(fireFreqDT)
  amc::.gc()

  ## add noFires to cohortData - no year info, because its the total across the simulation
  cols <- c("scenario", "rep", "pixelIndex", "noFires")
  allPixelCohortData <- unique(allPixelBurnData[, ..cols])[allPixelCohortData,
                                                           on = .(scenario, rep, pixelIndex)]

  ## checks
  if (mod$doAssertion) {
    cols <- c("scenario", "rep", "pixelIndex")
    test1 <- allPixelBurnData[allPixelCohortData[is.na(noFires), ..cols], on = cols, nomatch = 0]
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
  allPixelCohortData[, firePresAbs := as.integer(any(noFires > 0)), by = .(scenario, rep, pixelIndex)]
  amc::.gc()

  ## export to sim
  sim$allPixelCohortData <- allPixelCohortData
  sim$allPixelBurnData <- allPixelBurnData

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

addVegTypesCNEvent <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
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
    reps[reps < 10] <- paste0("0", reps[reps < 10])
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

  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

### add additional events as needed by copy/pasting from above
