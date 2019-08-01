# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects, use sim$xxx, and are thus globally available
# to all modules. Functions can be used without sim$ as they are namespaced, like functions
# in R packages. If exact location is required, functions will be: sim$<moduleName>$FunctionName
defineModule(sim, list(
  name = "Biomass_fuels",
  description = "SpaDES version of the LANDIS-II Dynamic Biomass Fuels Extention v2.2 - 15 Jun 2017", #"insert module description here",
  keywords = c("fire fuels", "fuel type", "LANDIS", "LandR"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca",
                   role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_fuels.Rmd"),
  reqdPkgs = list("data.table", "dplyr",
                  "PredictiveEcology/SpaDES.core@development",
                  "PredictiveEcology/SpaDES.tools@development",
                  "PredictiveEcology/reproducible@development",
                  "LandR"),
  parameters = rbind(
    defineParameter(name = "hardwoodMax", class = "numeric", default = 15L,
                    desc = "Threshold of percent biomass below which fuel types are considered conifer or mixed.
                    Defaults to 15, as in LANDIS example file"),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("successionTimestep", "numeric", 10L, NA, NA, "defines the simulation time step, default is 10 years")
    ),
  inputObjects = bind_rows(
    expectsInput(objectName = "sppEquiv", objectClass = "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput(objectName = "cohortData", objectClass = "data.table",
                 desc = "age cohort-biomass table hooked to pixel group map by pixelGroupIndex at
                 succession time step", sourceURL = NA),
    expectsInput(objectName = "species", objectClass = "data.table",
                 desc = "Species table produced by LandR-Biomass that has species traits such as longevity...",
                 sourceURL = NA),
    expectsInput(objectName = "speciesTable", objectClass = "data.table",
                 desc = "species attributes table, default is from Dominic and Yan's project",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput(objectName = "sppMultipliers", objectClass = "data.table",
                 desc = "Table of species biomass coefficient weights.
                 Recommended to be close to 1.0 for all species (see LANDIS-II Dynamic Fire System Extension (v2.1) User Guide).
                 Default values adapted from https://raw.githubusercontent.com/CeresBarros/Extension-Dynamic-Biomass-Fuels/master/testings/version-tests/v6.0-2.0/dynamic-biomass-fuels.txt"),
    expectsInput(objectName = "FuelTypes", objectClass = "data.table",
                 desc = "Table of Fuel Type parameters, with  base fuel type, species (in LANDIS code), their - or + contribution ('negSwitch'),
                 min and max age for each species (see LANDIS-II Dynamic Fire System Extension (v2.1) User Guide).
                 Default values adapted from https://raw.githubusercontent.com/CeresBarros/Extension-Dynamic-Biomass-Fuels/master/testings/version-tests/v6.0-2.0/dynamic-biomass-fuels.txt"),
    expectsInput(objectName = "fTypeEcore", objectClass = "data.table",
                 desc = "Table of Fuel Types per Ecoregion (optional, see LANDIS-II Dynamic Fire System Extension (v2.1) User Guide).
                 Default values adapted from https://raw.githubusercontent.com/CeresBarros/Extension-Dynamic-Biomass-Fuels/master/testings/version-tests/v6.0-2.0/dynamic-biomass-fuels.txt")
    ),
  outputObjects = bind_rows(
    createsOutput(objectName = "pixelFuelTypes", objectClass = "data.table",
                  desc = "Fuel types per pixel group, calculated from cohort biomasses")
  )
    ))

doEvent.Biomass_fuels = function(sim, eventTime, eventType) {
  switch(
    eventType,
    init = {
      # do stuff for this event
      sim <- fuelsInit(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, start(sim) + P(sim)$successionTimestep,
                           "Biomass_fuels", "doPrepareInputTables", eventPriority = 1)
      sim <- scheduleEvent(sim, start(sim) + P(sim)$successionTimestep,
                           "Biomass_fuels", "doFuelTypes", eventPriority = 1.5)
    },
    doPrepareInputTables = {
      # do stuff for this event
      sim <- prepareInputTables(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$successionTimestep,
                           "Biomass_fuels", "doPrepareInputTables", eventPriority = 1)
    },
    doFuelTypes = {
      # do stuff for this event
      sim <- calcFuelTypes(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, time(sim) + P(sim)$successionTimestep,
                           "Biomass_fuels", "doFuelTypes",  eventPriority = 1.5)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}

### template initialization
fuelsInit <- function(sim) {
  ## no module initalisation necessary
  return(invisible(sim))
}

prepareInputTables <- function(sim) {
  ## SPECIES COEFFICIENTS --------------------------------
  sppMultipliers <- copy(sim$sppMultipliers)

  ## check whether spp names match LANDIS or LandR format
  ## before merging species names and codes
  if (time(sim) == (start(sim) + sim@params$LBMR$successionTimestep)) {
    sppMultipliers <- prepSppMultipliers(sppMultipliers,
                                         sppEquiv = sim$sppEquiv,
                                         sppEquivCol = P(sim)$sppEquivCol)
  }

  ## FUEL TYPES TABLE ---------------------------------------------
  FuelTypes <- copy(sim$FuelTypes)

  ## check whether spp names match LANDIS or LandR format
  if (time(sim) == (start(sim) + sim@params$LBMR$successionTimestep)) {
    FuelTypes <- prepFuelTypes(FuelTypes,
                               sppEquiv = sim$sppEquiv,
                               sppEquivCol = P(sim)$sppEquivCol)

  }

  ## FUEL TYPES AND ECOREGIONS TABLE ----------------------
  ## assign NAs in fuel types with no ecoregion
  fTypeEcoreg <- copy(sim$fTypeEcoreg)
  fTypeEcoreg[FuelTypes[, .(FuelType)], on = "FuelType", nomatch = NA]
  fTypeEcoreg <- fTypeEcoreg[!duplicated(fTypeEcoreg)]


  ## export to sim
  sim$sppMultipliers <- sppMultipliers
  sim$FuelTypes <- FuelTypes
  sim$fTypeEcoreg <- fTypeEcoreg

  return(invisible(sim))
}

calcFuelTypes <- function(sim) {
browser()
    ## PIXEL FUEL TYPES TABLE ------------------------
  ## create pixelFuelTypes table from cohorData
  ## subset cohort data and non-na fuel types
  pixelFuelTypes <- copy(sim$cohortData[, pixelGroup:B])
  tempFT <- na.omit(copy(sim$FuelTypes[, -c("FuelTypeDesc"), with = FALSE]))  ## keep only complete lines with spp codes

  ## merge the two tables
  pixelFuelTypes <- pixelFuelTypes[tempFT, on = .(speciesCode), allow.cartesian = TRUE, nomatch = 0] %>%
    .[order(pixelGroup)]
  ## add sppMultipliers
  pixelFuelTypes <- pixelFuelTypes[sim$sppMultipliers[,.(speciesCode, Coefficient)],
                                   on = .(speciesCode), nomatch =0] %>%
    .[order(pixelGroup)]

  ## add fuel type ecoregions
  pixelFuelTypes <- pixelFuelTypes[sim$fTypeEcoreg, on = .(FuelType), nomatch = 0] %>%
    .[, ftEcoregion := as.numeric(Ecoregions)] %>%
    .[, Ecoregions := NULL]


  ## CHECK FUEL TYPES AND ECOREGIONS ------------------------------
  ## if fuel types are ecorgeion-specific, remove fuel types
  ## that are in wrong ecoregion
  if (any(!is.na(pixelFuelTypes$ftEcoregion))) {
    subsetDT <- pixelFuelTypes[!is.na(ftEcoregion)] %>%
      .[ftEcoregion == ecoregionGroup]
    subsetDT <- rbind(pixelFuelTypes[is.na(ftEcoregion)], subsetDT)
    pixelFuelTypes <- subsetDT
  }

  ## CALCULATE SPP VALUES FOR EACH FUEL TYPE IN EACH PIXEL ------------------------------
  ## calculate total biomass per pixelGroup, ecoregion and fuel type
  ## only species contributing to a given fuel type and with appropriate age are considered
  ## species biomass is weighted by the coeff and becomes
  ## negative if the species has a negative contribution (negSwitch) to that fuel type
  cols <- c("pixelGroup", "ecoregionGroup", "FuelType")
  pixelFuelTypes <- pixelFuelTypes[age >= minAge & age <= maxAge,
                                   forTypValue := sum(B*Coefficient*negSwitch),
                                   by = cols]
  ## remove lines that have no biomass for a fuel
  pixelFuelTypes <- pixelFuelTypes[!is.na(forTypValue)]

  ## ASSESS DOMINANT FUEL TYPE ----------------------------------------
  ## get max spp value (total biomass) in each pixelGroup and
  ## attribute corresponding fuel type in function of conifer/deciduous biomass
  cols <- c("pixelGroup", "ecoregionGroup")

  pixelFuelTypes[, maxValue := max(forTypValue, na.rm = TRUE),
                 by = cols]

  pixelFuelTypes[, finalBaseFuel := as.character()]
  pixelFuelTypes[, c("sumConifer", "sumDecid",
                     "coniferDom", "hardwoodDom",
                     "finalFuelType") := as.integer()]

  ## TODO: this takes ages and needs to be optimized for speed.
  pixelFuelTypes[, c("sumConifer", "sumDecid",
                     "coniferDom", "hardwoodDom",
                     "finalBaseFuel", "finalFuelType") := calcFinalFuels(BaseFuel, FuelType,
                                                                         forTypValue, maxValue,
                                                                         P(sim)$hardwoodMax),
                 by = cols]
  ## remove unnecessary columns and export to sim
  pixelFuelTypes <- cbind(pixelFuelTypes[, .(pixelGroup, ecoregionGroup)],
                          pixelFuelTypes[, finalBaseFuel:finalFuelType]) %>%
    .[!duplicated(.)]
  sim$pixelFuelTypes <- copy(pixelFuelTypes)

  return(invisible(sim))
}

.inputObjects <- function(sim) {
  dPath <- dataPath(sim)
  ## Get LANDIS example parameters -----------------------
  ## to use if others haven't been supplied in <module>/inputs
  if (any(!suppliedElsewhere("sppMultipliers", sim),
          !suppliedElsewhere("FuelTypes", sim),
          !suppliedElsewhere("fTypeEcoreg", sim)) &
      any(!file.exists(file.path(dataPath(sim), "sppMultipliers.csv")),
          !file.exists(file.path(dataPath(sim), "FuelTypes.csv")),
          !file.exists(file.path(dataPath(sim), "fTypesEcoregions.csv")))) {
    maxcol <- 21 #max(count.fields(file.path(getPaths()$dataPath, "dynamic-biomass-fuels.txt"), sep = ""))
    dynamicBiomassFuels <- Cache(prepInputs,targetFile = "dynamic-biomass-fuels.txt",
                                 url = paste0("https://raw.githubusercontent.com/CeresBarros/",
                                              "Extension-Dynamic-Biomass-Fuels/master/testings/",
                                              "version-tests/v6.0-2.0/dynamic-biomass-fuels.txt"),
                                 destinationPath = dataPath(sim),
                                 fun = "utils::read.table",
                                 fill = TRUE, row.names = NULL,
                                 sep = "",
                                 header = FALSE,
                                 blank.lines.skip = TRUE,
                                 col.names = c(paste("col",1:maxcol, sep = "")),
                                 stringsAsFactors = FALSE)
    dynamicBiomassFuels <- data.table(dynamicBiomassFuels[, 1:15])
    dynamicBiomassFuels <- dynamicBiomassFuels[!(col1 %in% c("LandisData", "Timestep", "MapFileNames",
                                                             "PctConiferFileName", "PctDeadFirFileName"))]
    dynamicBiomassFuels <- dynamicBiomassFuels[!(col2 %in% c("The", "Users", "Optional"))]
    dynamicBiomassFuels <- dynamicBiomassFuels[!(col1 == ">>" & grepl("---", col2))]
    dynamicBiomassFuels[col1 == ">>"] <- data.table(dynamicBiomassFuels[col1 == ">>", col2:col15], col15 = "NA")
    dynamicBiomassFuels[,col15:=NULL]
    dynamicBiomassFuels[col1 == "Fuel" & col2 == "Type"] <- data.table(dynamicBiomassFuels[col1 == "Fuel" & col2 == "Type", col2:col14],
                                                                       col14 = "NA")
    dynamicBiomassFuels[col1 == "Type", col1 := "FuelType"]
  }

  ## SPECIES COEFFICIENTS ---------------------------------
  if (!suppliedElsewhere("sppMultipliers", sim)) {
    if (file.exists(file.path(dataPath(sim), "sppMultipliers.csv"))) {
      sim$sppMultipliers <- prepInputs(targetFile = "sppMultipliers.csv",
                                       destinationPath = dataPath(sim),
                                       fun = "read.csv",
                                       header = TRUE) %>%
        data.table()
    } else {
      message(paste0("Can't find sppMultipliers.csv in ", dataPath(sim),
                     ".\nUsing LANDIS example file"))

      sppMultipliers <- dynamicBiomassFuels[(which(col1=="Fuel") + 1) : (which(col1 == "HardwoodMaximum") - 1),
                                            col1:col2]

      names(sppMultipliers) <- as.character(sppMultipliers[1,])
      sppMultipliers <- sppMultipliers[-1]

      sim$sppMultipliers <- copy(sppMultipliers)
    }
  }

  ## FUEL TYPES TABLE -------------------------------------
  if (!suppliedElsewhere("FuelTypes", sim)) {
    if (file.exists(file.path(dataPath(sim), "FuelTypes.csv"))) {
      sim$FuelTypes <- prepInputs(targetFile = "FuelTypes.csv",
                                  destinationPath = dataPath(sim),
                                  fun = "utils::read.csv",
                                  header = TRUE) %>%
        data.table

    } else {
      message(paste0("Can't find FuelTypes.csv in ", dataPath(sim),
                     ".\nUsing LANDIS example file"))

      FuelTypes <- dynamicBiomassFuels[(which(col1=="FuelTypes") + 1) : (which(col1 == ">>EcoregionsTable") - 1),
                                       col1:col14]
      ## rename columns
      FuelTypes[1, `:=`(col3 = "minAge",
                        col4 = "NA",
                        col5 = "maxAge",
                        col6 = "Species")]
      FuelTypes[, col4 := NULL]
      names(FuelTypes) <- c(as.character(FuelTypes[1, 1:5, with = FALSE]), paste0("Species", 2:9))
      FuelTypes <- FuelTypes[-1]

      ## melt table to get all species in one column
      FuelTypes <- melt(FuelTypes, id.vars = 1:4, variable.name = "Species") %>%
        .[value != ""] %>%
        .[, Species := value] %>%
        .[, value := NULL]

      ## create a column with negative switch
      ## species with a negative switch negatively contribute to a fuel type
      FuelTypes[grepl("-", Species), negSwitch := -1L]
      FuelTypes[is.na(negSwitch), negSwitch := 1]
      FuelTypes[, Species := sub("-", "", Species)]

      ## export table to sim
      sim$FuelTypes <- copy(FuelTypes)
    }
  }

  ## FUEL TYPES AND ECOREGIONS TABLE ----------------------
  if (!suppliedElsewhere("fTypeEcoreg", sim)) {
    if (file.exists(file.path(dataPath(sim), "fTypesEcoregions.csv"))) {
      sim$fTypeEcoreg <- prepInputs(targetFile = "fTypesEcoregions.csv",
                                    destinationPath = dataPath(sim),
                                    fun = "utils::read.csv",
                                    header = TRUE) %>%
        data.table

    } else {
      message(paste0("Can't find fTypesEcoregions.csv in ", dataPath(sim),
                     ".\nUsing LANDIS example file"))
      ## in the LANDIS test exeamples this feature is turned off and the
      ## ecoregions do not correspond to those in the succession module examples.
      ## assigning NAs for ecoregion sin all fuel types to simualte this behaviour
      fTypeEcoreg <- data.table(FuelType = sort(unique(sim$FuelTypes$FuelType)), Ecoregions = NA)
      sim$fTypeEcoreg <- fTypeEcoreg
    }
  }

  ## SPECIES EQUIVALENCY TABLE ---------------------------
  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)

    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                               newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide please provide sppColorVect")
  }

  ## SPECIES TRAITS TABLE ---------------------------------
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- prepInputs("speciesTraits.csv",
                                   destinationPath = dPath,
                                   url = extractURL("speciesTable"),
                                   fun = "utils::read.csv",
                                   header = TRUE, stringsAsFactors = FALSE) %>%
      data.table()
  }



  return(invisible(sim))
}
