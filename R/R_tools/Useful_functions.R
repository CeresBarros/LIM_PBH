## ------------------------------------------------------
## USEFUL FUNCTIONS
##
## Ceres: March 2020
## ------------------------------------------------------

## DEFINE FIRE EVENTS -----------------------
## wrapper around .calculateFireEvents

## Method from Andison (2012), defines a fire event composed of
## disturbed patches (severity/mortality >= 95%)
## island remnants (severity < 95%, surrounded by disturbed patches)
## matrix remnants (undisturbed, and partially surrouned by disturbed patches)

## sfObj is a Simple Features object from sf package
## fireNAMES should be a column/attribute with fire ID/names
## fireVARS can be NULL, or a vector of variable names/indices to retain
## crsProj is a string of the projection to use. If NULL (default) will use the same projection as sfObj, otherwise sfObj will be reprojected
## buff.dist if the buffer distance to define fire events
## PLOT, SAVE and overwrite determine if plotting, saving and overwriting should be done
## outputDIR and fileNAME define the directory and fileNAME to save the fire events shapefile (".shp" will be added to the name string),
## outputDIR will be created if non-existent

defineFireEvents <- function(sfObj, fireNAMES = NULL, fireVARS = NULL, crsProj = NULL, buff.dist = NULL, PLOT = TRUE, SAVE = TRUE,
                             outputDIR = NULL, fileNAME = NULL, overwrite = TRUE) {
  ## checks
  if(any(class(sfObj) != c("sf", "data.frame"))) stop ("sfObj must be an sf object")
  if(is.null(buff.dist)) stop("Define a buffer distance")
  if(is.null(fireNAMES)) stop("Provide the name of the fire ID variable")
  if(SAVE & is.null(outputDIR)) stop("SAVE is TRUE, but output folder is not defined")
  if(SAVE & is.null(fileNAME)) stop("SAVE is TRUE, but file name prefix is not defined")

  ## DEFINE PROJECTION AND RE-PROJECT IF NEED BE
  crsProj <- if (is.null(crsProj)) crs(sfObj) else CRS(crsProj)

  if (!compareCRS(crsProj, crs(sfObj))) {
    warning("Reprojecting sfObj to selected projection")
    sfObj <- st_transform(sfObj, crs = st_crs(crsProj))
  }

  ## Get vector of fire names
  fire.ls <- unique(eval(parse(text = paste0("sfObj$", fireNAMES))))

  ## CHECK FIRE VARIABLES
  if(is.numeric(fireVARS)) fireVARS <- names(sfObj)[fireVARS]

  if(!is.null(fireVARS)) cat(paste0("Using the following fire variables: ", paste0(fireVARS, collapse =", ")), "\n")

  ## CALCULATE FIRE EVENTS
  fireEvent.ls <- lapply(fire.ls, FUN = .calculateFireEvents,
                         sfObj = sfObj, fireNAMES = fireNAMES, fireVARS = fireVARS,
                         crsProj = crsProj, buff.dist = buff.dist)

  fireEvent.all <- do.call(rbind, fireEvent.ls)

  if(PLOT) plot(fireEvent.all, key.pos = 1)

  ## SAVE AS SHAPEFILE
  if(SAVE) {
    ## clean ws before saving
    rm(fireEvent.ls, sfObj); gc(reset = TRUE)

    if(!dir.exists(outputDIR)) dir.create(outputDIR, recursive = TRUE)
    st_write(fireEvent.all, file.path(outputDIR, paste0(fileNAME, ".shp")), delete_layer = overwrite)
  }

  return(fireEvent.all)
}

## CALCULATE FIRE EVENTS -----------------------
## Method from Andison (2012), defines a fire event composed of
## disturbed patches (severity/mortality >= 95%)
## island remnants (severity < 95%, surrounded by disturbed patches)
## matrix remnants (undisturbed, and partially surrouned by disturbed patches)

## sfObj is a Simple Features object from sf package
## fireNAMES should be a column/attribute with fire ID/names
## fireVARS can be NULL, or a vector of variable names/indices to retain
## crsProj is a string of the projection to use. If NULL (default) will use the same projection as sfObj, otherwise sfObj will be reprojected
## buff.dist if the buffer distance to define fire events

.calculateFireEvents <- function(fire, sfObj, fireNAMES, fireVARS,
                                 crsProj, buff.dist) {
  print(as.character(fire))

  firePolys <- eval(parse(text = paste0("sfObj$", fireNAMES))) == fire

  if (is.null(fireVARS)) {
    sf.fire <- sfObj[firePolys, c(fireNAMES)]
  } else sf.fire <- sfObj[firePolys, c(fireNAMES, fireVARS)]

  ## CALCULATE FIRE AND EVENT PERIMETERS - buffer out and then in.
  firePerim <- st_union(sf.fire$geometry)
  outerBuff <- st_buffer(firePerim, dist = buff.dist)
  eventPerim <- st_buffer(outerBuff, dist = -buff.dist)

  ## REMOVE INNER MATRIX HOLES FROM EVENT AND ORIGINAL FIRE PERIMETER
  ## (i.e. unburnt patches surrounded by fire)
  if(class(eventPerim)[1] == "sfc_MULTIPOLYGON") {
    noHolesEventPerim <- st_sfc(st_multipolygon(lapply(eventPerim[[1]], function(x) x[1])), crs = st_crs(crsProj))
    noHolesEventPerim <- st_union(noHolesEventPerim)  ## union to account for disturbed patches inside inner matrix, nested within larger disturbed patches
  } else {
    noHolesEventPerim <- st_sfc(st_multipolygon(lapply(eventPerim[1], function(x) x[1])), crs = st_crs(crsProj))
    noHolesEventPerim <- st_union(noHolesEventPerim) ## union to account for disturbed patches inside inner matrix, nested within larger disturbed patches
  }

  if(class(firePerim)[1] == "sfc_MULTIPOLYGON") {
    noHolesFirePerim <- st_sfc(st_multipolygon(lapply(firePerim[[1]], function(x) x[1])), crs = st_crs(crsProj))
    noHolesFirePerim <- st_union(noHolesFirePerim)
  } else {
    noHolesFirePerim <- st_sfc(st_multipolygon(lapply(firePerim[1], function(x) x[1])), crs = st_crs(crsProj))
    noHolesFirePerim <- st_union(noHolesFirePerim)
  }

  ## EXTRACT INNER MATRIX HOLES PRODUCED BY BUFFERING - these will determine what the inner matrix remnants are
  bufferedHoles <- st_difference(noHolesEventPerim, eventPerim)

  ## EXTRACT ALL REMNANTS PRODUCED BY BUFFERING (except holes)
  allResiduals <- st_difference(eventPerim, firePerim)

  ## CALCULATE INTERSECTIONS (TOUCH) BETWEEN HOLES AND RESIDUALS
  if(length(bufferedHoles) > 0) {
    if(class(allResiduals)[1] == "sfc_MULTIPOLYGON") {
      remnHoleInters  <- sapply(allResiduals[[1]], FUN = function(sfgpoly) {
        sfcpoly <- st_sfc(list = st_polygon(sfgpoly[1]), crs = st_crs(crsProj))
        st_intersects(sfcpoly, bufferedHoles, sparse = FALSE)
      })
    } else {
      remnHoleInters  <- sapply(allResiduals[1], FUN = function(sfgpoly) {
        sfcpoly <- st_sfc(list = st_polygon(sfgpoly[1]), crs = st_crs(crsProj))
        st_intersects(sfcpoly, bufferedHoles, sparse = FALSE)
      })
    }
  } else {
    remnHoleInters <- if(class(allResiduals)[1] == "sfc_MULTIPOLYGON") {
      rep(FALSE, length(allResiduals[[1]]))
    } else {
      rep(FALSE, length(allResiduals[1]))
    }
  }

  ## DEFINE INTERIOR RESIDUALS
  if(class(allResiduals)[1] == "sfc_MULTIPOLYGON") {
    inResids <- lapply(allResiduals[[1]][remnHoleInters], st_polygon) %>% st_sfc(., crs = st_crs(crsProj))
    inResids <- st_union(inResids)
  } else {
    inResids <- lapply(allResiduals[1][remnHoleInters], st_polygon) %>% st_sfc(., crs = st_crs(crsProj))
    inResids <- st_union(inResids)
  }

  ## DEFINE OUTER MATRIX REMNANTS
  outResids <- st_difference(allResiduals, inResids)

  ## MERGE HOLES AND INNER MATRIX
  inResids2 <- c(bufferedHoles, inResids)

  ## CONVERT TO MULTIPOLYGONS AND SIMPLE FEATURE COLLECTION, ADDING FIRE NAME
  firePerim <- st_sfc(firePerim)
  firePerim <- st_sf(geometry = firePerim)   ## first "combine" list of polygons into a multipolygon, which is then converted to a Simple Features object
  firePerim$PatchType <- "disturbedPatch"

  eventPerim <- st_sfc(eventPerim, check_ring_dir = TRUE)
  eventPerim <- st_sf(geometry = eventPerim)
  eventPerim$PatchType <- "eventPerim"

  outResids <- st_sfc(outResids)
  outResids <- st_sf(geometry = outResids)
  outResids$PatchType <- "outResids"

  inResids2 <- st_sfc(inResids2)
  inResids2 <- st_sf(geometry = inResids2)
  if (length(inResids2$geometry))
    inResids2$PatchType <- "inResids"

  ## add fire details
  if(!is.null(fireVARS)) {
    firePerim <- st_join(firePerim, sf.fire, left = FALSE)     ## intersection because firePerim is entirely within sf.fire
    eventPerim <- st_join(eventPerim, sf.fire, left = FALSE)   ## using inner join (see ?inner_join)
    temp.df <- firePerim[, !names(firePerim) %in% names(outResids), drop = TRUE]

    if("SEV_CLASS" %in% fireVARS) temp.df$SEV_CLAS <- NA

    temp.df <- temp.df[!duplicated(temp.df),]

    outResids <- merge(outResids, temp.df)
    inResids2 <- merge(inResids2, temp.df)

    ## COMBINE INTO ONE OBJECT
    fireEvent <- rbind(firePerim, eventPerim, outResids, inResids2)
  } else {
    ## COMBINE INTO ONE OBJECT
    fireEvent <- rbind(firePerim, eventPerim, outResids, inResids2)

    ## ADD FIRE NAME
    eval(parse(text = paste0("fireEvent$", fireNAMES, "<- fire")))
  }

  return(fireEvent)
}

## CALCULATE NEIGHBOURHOOD SEVERITY -----------------------
## Calculates the average severity of neighbours according
## to a given distance or set of distances

## dist is a vector of distances in meters - note that sevPoints much be in a meter-based projection
## sevPoints is an sf object of points and severity
## sevColID is the columns name of the severity column

calculateNgbSevWrapper <- function(dists, sevPoints, sevColID, parallel = TRUE) {
  if (length(dists) > 1) {
    if (parallel) {
      message("Starting parallelization...")
      require(future.apply)
      plan(multisession(gc = TRUE))
      ngbSEVList <- future_lapply(dists, FUN = .calculateNgbSev,
                                  sevPoints = sevPoints, sevColID = sevColID)

    } else {
      ngbSEVList <- lapply(dists, FUN = .calculateNgbSev,
                           sevPoints = sevPoints, sevColID = sevColID)
    }
    ngbSEVDT <- Reduce(function(x, y) merge(x, y, by = "pixID", all = TRUE),
                       ngbSEVList)
  } else
    ngbSEVDT <- .calculateNgbSev(dists, sevPoints, sevColID)

  return(ngbSEVDT)
}

.calculateNgbSev <- function(dist, sevPoints, sevColID) {
  if (sum(names(sevPoints) %in% sevColID) > 1)
    stop("Several column names match 'sevColID")
  if (!sum(names(sevPoints) %in% sevColID))
    stop("No column names match 'sevColID")

  ## change col name for generality
  names(sevPoints) <- sub(sevColID, "sev", names(sevPoints))

  ## draw buffers
  message(paste0("Drawing ", dist, "m buffers and detecting neighbours"))
  bufferSf <- st_buffer(sevPoints, dist = dist) ## keep all columns so that the join identifies .x and .y columns

  ## join to find with pixels are within another's buffer
  ## st_touches avoids "self" joins, so pixels are not joined with their own buffer
  pointsWithinBuffer <- st_join(bufferSf, sevPoints, join = st_touches)
  names(pointsWithinBuffer) <- sub("\\.x", "buffer", names(pointsWithinBuffer))
  names(pointsWithinBuffer) <- sub("\\.y", "points", names(pointsWithinBuffer))

  pointsWithinBufferDT <- data.table(st_drop_geometry(pointsWithinBuffer))
  pointsWithinBufferDT[, sevbuffer := NULL] ## keep track of neighbour sev only

  setnames(pointsWithinBufferDT, old = c("pixIDbuffer", "pixIDpoints", "sevpoints"),
           new = c("pixID", "pixIDneigh", "sevngb"))

  message(paste0("Calculating average severity across neighbours"))
  ngbhoodSEV <- pointsWithinBufferDT[, list(sevngbhood = mean(sevngb, na.rm = TRUE)),
                                     by = pixID]
  setnames(ngbhoodSEV, old = "sevngbhood",
           new = paste0("meanngb", sevColID, "_", dist, "m"))
  message(paste0("Done!"))
  return(ngbhoodSEV)
}


## CALCULATE NEIGHBOURHOOD NO. BURNT PIXELS -----------------------
## Calculates the proportion of burnt neighbours according
## to a given distance or set of distances. Points with 0 severity are assumed to not be burnt.

## dist is a vector of distances in meters - note that sevPoints much be in a meter-based projection
## sevPoints is an sf object of points and severity
## sevColID is the columns name of the severity column
## cores is the number of cores to use for parallelisation

calculateNgbBurnsWrapper <- function(dists, sevPoints, sevColID, fireColID,
                                     resolution = resolution, parallel = TRUE, cores = NULL) {
  ## make a list of combinations of fire ID and buffer distance
  fireBufferCombos <- expand.grid(unique(sevPoints[[fireColID]]), dists)
  names(fireBufferCombos) <- c("fireID", "dists")

  if (nrow(fireBufferCombos) > 1) {
    if (parallel) {
      message("Starting parallelization...")
      require(future.apply)
      if (is.null(cores)) {
        cores <- availableCores()
      }
      plan(multisession(gc = TRUE), workers = cores)
      ngbSEVList <- future_mapply(FUN = .calculateNgbBurns, dist = fireBufferCombos$dists,
                                  fireID = fireBufferCombos$fireID,
                                  MoreArgs = list(sevPoints = sevPoints, sevColID = sevColID,
                                                  fireColID = fireColID, resolution = resolution),
                                  SIMPLIFY = FALSE)

    } else {
      ngbSEVList <- Cache(Map,
                          dist = fireBufferCombos$dists,
                          fireID = fireBufferCombos$fireID,
                          MoreArgs = list(sevPoints = sevPoints, sevColID = sevColID,
                                          fireColID = fireColID, resolution = resolution),
                          .calculateNgbBurns)
    }
    ngbSEVDT <- rbindlist(ngbSEVList)
  } else
    ngbSEVDT <- .calculateNgbBurns(fireBufferCombos$dists, fireBufferCombos$fireID,
                                   sevPoints, sevColID, fireColID, resolution)

  return(ngbSEVDT)
}

.calculateNgbBurns <- function(dist, fireID, sevPoints, sevColID, fireColID, resolution) {
  if (sum(names(sevPoints) %in% sevColID) > 1)
    stop("Several column names match 'sevColID")
  if (!sum(names(sevPoints) %in% sevColID))
    stop("No column names match 'sevColID")

  ## change col name for generality
  names(sevPoints) <- sub(sevColID, "sev", names(sevPoints))

  ## do one fire at a time
  ## make raster of severity and raster of pixel IDs
  i <- which(sevPoints[[fireColID]] == fireID)
  fireRas <- raster(as_Spatial(sevPoints[i,]), resolution = resolution, crs = crs(sevPoints))
  fireRas[] <- NA
  fireRas <- rasterize(as_Spatial(sevPoints[i,]), fireRas, field = "sev")
  fireRasIDs <- rasterize(as_Spatial(sevPoints[i,]), fireRas, field = "pixID")

  ## draw buffers
  message(paste0("Drawing ", dist, "m buffers and counting burnt neighbours for ", fireID, " fire"))
  w <- focalWeight(fireRas, d = dist*3, type = "circle")
  w[w > 0] <- 1
  w[w == 0] <- NA
  w[ceiling(nrow(w)/2), ceiling(ncol(w)/2)] <- 0  ## exclude focal cell

  ## calculate prop burnt neighbours per pixel
  message(paste0("Calculating number of burnt neighbours"))
  ngbhoodBurns <- focal(fireRas, w = w, pad = TRUE,
                        fun = function(x) {sum(x > 0, na.rm = TRUE)/sum(!is.na(x))})
  ngbhoodBurns <- mask(ngbhoodBurns, fireRas)  ## need to remove 0s beyond fire perimeter

  ## make DT
  ngbhoodBurnsDT <- data.table(pixID = getValues(fireRasIDs),
                               ngbPropBurns = getValues(ngbhoodBurns),
                               fire = fireID,
                               bufferSize = dist)
  ngbhoodBurnsDT <- na.omit(ngbhoodBurnsDT)
  setnames(ngbhoodBurnsDT, old = "fire", new = fireColID)

  ## checks
  if (length(unique(ngbhoodBurnsDT$pixID)) != length(i)) {
    noMissing <- sum(!sevPoints$pixID[i] %in% ngbhoodBurnsDT$pixID)
    warning(paste(noMissing, "points in fire perimeter were not converted to raster pixels.",
                  "\nThis is probably due to more than one point falling in the same cell"))
  }
  message(paste0("Done!"))
  return(ngbhoodBurnsDT)
}



## WRAPPER FUNCTION TO ESTIMATE HYPERVOLUME BANDWIDTHS  -----------------------
##
## allData is a data.table with the data for both hypervolumes and an ID column (HVidvar) fo
##   that identifies the data for each hypervolume
## ... further arguments passed to ToolsCB::HVordination

estimateBW_wrapper <- function(allData, HVidvar, ...) {
  HVnames <- unique(allData[, HVidvar])

  ## need to re-do categorical variables so that levels correspond to unique values
  ordi.list <- HVordination(datatable = allData, HVidvar = HVidvar, ...)

  HVpoints <- ordi.list[[1]]
  noAxes <- ordi.list[[2]]

  HV1rows <- allData[, HVidvar] %in% HVnames[1]
  HV2rows <- allData[, HVidvar] %in% HVnames[2]

  temp <- data.frame(SilvBW_HV1 = hypervolume::estimate_bandwidth(HVpoints[HV1rows, 1:noAxes]),
                     SilvBW_HV2 = hypervolume::estimate_bandwidth(HVpoints[HV2rows, 1:noAxes]),
                     stdev_HV1 = apply(HVpoints[HV1rows, 1:noAxes], 2, sd),
                     stdev_HV2 = apply(HVpoints[HV2rows, 1:noAxes], 2, sd),
                     PC = c(1:noAxes))
  temp$HVpair = paste0(HVnames[1], "_", HVnames[2])
  return(temp)
}

## GET PCA LOADINGS -----------------------
## code from \code{biplot.prcomp} and \code{biplot.default}
## x is \code{prcomp} object
getLoadings4Plot <- function(x, choices = c(1,2), scale = 1, pc.biplot = FALSE,
                             xlim, ylim, expand = 1) {
  if (!length(scores <- x$x[,choices])) {
    stop(gettextf("object '%s' has no scores", deparse1(substitute(x))),
         domain = NA)
  }
  if (is.complex(scores)) {
    stop("biplots are not defined for complex PCA")
  }
  lam <- x$sdev[choices]
  n <- NROW(scores)
  lam <- lam * sqrt(n)
  if (scale < 0 || scale > 1) {
    warning("'scale' is outside [0, 1]")
  }
  if (scale != 0) {
    lam <- lam^scale
  } else {
    lam <- 1
  }
  if (pc.biplot) {
    lam <- lam/sqrt(n)
  }
  loadings <- t(t(x$rotation[,choices]) * lam)
  scores <- t(t(scores[, choices])/lam)

  ## rescale the loadings to plot within the PCA scores plot
  ## from biplot.default
  unsigned.range <- function(x) c(-abs(min(x, na.rm = TRUE)),
                                  abs(max(x, na.rm = TRUE)))
  rangx1 <- unsigned.range(scores[, 1L])
  rangx2 <- unsigned.range(scores[, 2L])
  rangy1 <- unsigned.range(loadings[, 1L])
  rangy2 <- unsigned.range(loadings[, 2L])
  if (missing(xlim) && missing(ylim)) {
    xlim <- ylim <- rangx1 <- rangx2 <- range(rangx1, rangx2)
  } else {
    if (missing(xlim)) {
      xlim <- rangx1
    } else {
      if (missing(ylim)) {
        ylim <- rangx2
      }
    }
  }

  ratio <- max(rangy1/rangx1, rangy2/rangx2)/expand
  loadings <- loadings/ratio
  return(loadings)
}



## SAMPLE SIMULATION YEARS -----------------------
#' SAMPLE SIMULATION YEARS
#' samples 5 years (per rep) at regular intervals (every century) within the last 500 years of sampling
#'
#' @param yearRepTable a table of years and repetitions (unique combos will be extracted)
#' @param .seed a numeric passed to `set.seed`. If NA, seed won't be set.
#'
#' @return a table with years to sample per rep


sample5SimYears <- function(yearRepTable, .seed = 123) {
  yearSamples <- setkeyv(unique(yearRepTable), c("rep", "year"))
  yearSamples[, group := cut(year, breaks = 5, right = FALSE, labels = FALSE), by = rep]

  if (!is.na(.seed)) {
    initialRandomSeed <- .Random.seed
    set.seed(.seed)
  }

  yearSamples[, year2 := sample(year, 1), by = .(rep, group)]
  needsNewSample <- yearSamples[, length(unique(year2)) < 5, by = group]
  while(any(needsNewSample$V1)) {
    yearSamples[group %in% needsNewSample[which(V1), group], year2 := sample(year, 1), by = .(rep, group)]
    needsNewSample <- yearSamples[, length(unique(year2)) < 5, by = group]
  }
  yearSamples <- unique(yearSamples[,.(year2, rep)])
  setnames(yearSamples, "year2", "year")

  if (exists("initialRandomSeed", inherits = FALSE)) {
    .Random.seed <- initialRandomSeed
  }
  return(yearSamples)
}


## Average AllPixelCDMntEnd columns across sampled years -----------------------
#' Average AllPixelCDMntEnd columns across sampled years
#'
#' Averages, takes the unique values, or the most frequent values
#' across years (per rep, pixel and species)
#' for all colums in allPixelCDMntEnd, except pixelGroup, which is
#' ignored and therefore excluded
#'
#' This function exists for caching purposes.

averagAllPixelCDMntEnd <- function(allPixelCDMntEnd) {
  allPixelCDMntEnd <- allPixelCDMntEnd[, list(
    vegType = as.integer(names(which.max(table(vegType)))),
    vegTypeCN = names(which.max(table(vegTypeCN))),
    noFires = as.integer(mean(noFires)),
    ecozoneCode = as.integer(unique(ecozoneCode)),
    ecozoneName = unique(ecozoneName),
    ecoregionGroup = unique(ecoregionGroup),
    age = as.integer(mean(age)),
    B = as.integer(mean(B)),
    mortality = as.integer(mean(mortality)),
    aNPPAct = as.integer(mean(aNPPAct)),
    firePresAbs = as.integer(names(which.max(table(firePresAbs))))
  )
  , by = .(scenario, rep, pixelIndex, speciesCode)]

  allPixelCDMntEnd
}


## XGBoost wrapper -----------------------
#' Wrapper for XGBoost
#'
#' Tune, fits and tests XGBoost models with
#' k-fold cross-validation.
#'
#' @param dat a `data.table` containing predictors and response variable.
#' @param dig a digest passed to `Cache(..., .cacheExtra)` to bypass digesting
#'    `dat`. Often `dig` is a digest of `dat`.
#' @param nFolds number of folds for cross-validating the final model
#'    (i.e. the model using tuned parameters).
#' @param colnamesResp Name of column in `dat` to use as response variable.
#'    All other columns will be used as predictors
#' @param interaction_constraints passed to `xgboost::xgboost`.
#'    By default no interaction constraints.
#' @param SHAPthresh. Quantile threshold used for feature (i.e. variable) selection
#'    based on SHAP values. Features with SHAP values below the quantile threshold
#'    are excluded and the model re-run. A warning is issued if this resulted in poorer
#'    performace (based on AUC score), in which case one may consider relaxing (i.e. lowering)
#'    the threshold.
#' @param figDir if not `NULL`, diagnostic tuning plots will be saved to this directory.
#'
#' @return a list (one entry per fold) of lists with:
#'   * `$mod`: fitted model
#'   * `shap_values`: SHAP values for the fitted model
#'   * `shap_long`: long version of the SHAP values for the fitted model (used for plotting)
#' @importFrom caret createFolds
#' @importFrom purrr pmap
#' @importFrom pROC roc
#' @importFrom crayon cyan
#' @importFrom SHAPforxgboost shap.values
#' @importFrom reproducible Cache
runXGBOOST <- function(dat, dig, nFolds = 5, colnamesResp = "SEV_PROP",
                       interaction_constraints = NULL, SHAPthresh = 0,
                       figDir = NULL) {
  # dat <- dat[sample(NROW(dat), size = 1e4)]
  # tt <- table(dat$SEV_PROP)
browser()
  # Add dummy variables for factor columns -- i.e., the random effects
  if (all(sapply(dat, is.numeric)) %in% FALSE)
    dat <- model.matrix(~ . + 0, data = dat) |>
      Cache(omitArgs = c("object", "data", "x")
            , .cacheExtra = dig
            # , cacheid = "f5c41c45b1ca9cbc"
            ) # Creates dummy variables

  dat <- as.data.table(dat)

  colnamesPred <- setdiff(colnames(dat), colnamesResp) ## after model.matrix bcs colnames change

  ## Setup k-folds -----
  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  yearColname <- grep("year", tolower(colnames(dat)), value = TRUE)
  indexNames <- c("allData", "evalData")

  if (length(yearColname)) {
    crossValType <- "time-ordered"
    times <- unique(dat[[yearColname]])
    testLength <- 3
    initialWindow <- length(times) - testLength - nFolds + 1
    trainIndexK <- createTimeSlices(times, initialWindow = initialWindow, testLength, fixedWindow = FALSE)
    trainIndexK <- Map(tr = trainIndexK$train, te = trainIndexK$test, function(tr, te) {
      allData <- which(dat[[yearColname]] %in% times[c(tr, te)])
      evalData <- which(dat[[yearColname]] %in% times[te])
      list(allData, evalData) |> setNames(indexNames)
    })
  } else {
    crossValType <- "crossValidation"
    ## create folds and make a list with indices of full dataset and eahc fold
    trainIndexK <- createFolds(dat[[colnamesResp]], k = nFolds, list = TRUE, returnTrain = FALSE)
    trainIndexK <- Map(tr = trainIndexK, function(tr) {
      list(seq(NROW(dat)), tr) |> setNames(indexNames)
    })
  }

  ## sample columns after setting seed for caching (if different, then cache is triggered)
  colOrder <- setdiff(colnames(dat), c(yearColname))
  colOrder <- sample(colOrder)
  dat <- dat[, ..colOrder]
  dig <- .robustDigest(dat)

  ## Tune parameters on full data with caret first ----
  params <- .tunexgboost(dig,
                         dat[, .SD, .SDcols = c(colnamesPred, colnamesResp)],
                         colnamesResp = colnamesResp,
                         figDir) |>
    Cache(omitArgs = c("dat", "figDir")
          #, cacheId = "7b1de6940b00ff51"
    )

  ## subset predictor data
  datPreds <- dat[, ..colnamesPred]

  st <- system.time(
    mm <- pmap(
      list(dataFolds = trainIndexK, kFold = seq(nFolds)),
      function(dataFolds, kFold) {
        ## get row IDs for training data (allData) and testing data
        allDataIDs <- dataFolds[[indexNames[[1]]]]
        testIDs <- dataFolds[[indexNames[[2]]]]   ## eval data
        dig2 <- .robustDigest(dataFolds)

        # xgboost objects do not save with `qs` ... must be `rds`
        # opt <- options(reproducible.cacheSaveFormat = "rds")
        # on.exit(options(opt)) # redundant; but necessary if it fails during fit
        lowSHAPcols <- 1
        calcThresh <- TRUE
        # SHAPthresh <- 0.25  ## test
        cols2keep <- colnames(datPreds)

        modOut <- NULL

        while (length(lowSHAPcols)) {
          ## TODO: test: go back to previous model if AUC decreases after removing features
          browser()

          modOut2 <- xgboost(x = datPreds[allDataIDs],
                             , y = dat[[colnamesResp]][allDataIDs]
                             , interaction_constraints = interaction_constraints
                             # , objective = "reg:tweedie" ## no improvements
                             , nthread = 10
                             , eval_set = testIDs,
                             , monitor_training = TRUE
                             , eval_metric = c("auc", "rmse", "logloss")
                             , early_stopping_rounds = 100
                             , max_depth = params$max_depth   ## improved fit.
                             , nrounds = params$nrounds
                             , learning_rate = params$eta
                             , min_split_loss = params$gamma
                             , min_child_weight = params$min_child_weight
                             , colsample_bytree = params$colsample_bytree

          ) |>
            Cache(omitArgs = c("x", "y", "eval_set"),
                  .functionName = .functionNameHelper("xgboost", kFold),
                  .cacheExtra = c(dig, dig2, cols2keep),
                  showSimilar = TRUE,
                  cacheSaveFormat = "rds")

          if (is.null(modOut)) {
            modOut <- modOut2
            AUCout <- tail(attr(modOut, "evaluation_log"), 1)$train_auc
          }

          ## get last AUC
          AUCout2 <- tail(attr(modOut2, "evaluation_log"), 1)$train_auc
          if (AUCout2 < AUCout) {
            message("AUC decreased after removing features.\n",
                    "  The previous model will be retained, instead")
          } else {
            modOut <- modOut2
          }
          AUCout <- AUCout2

          ## calculate predictions and residuals
          valData <- datPreds[testIDs,]
          valData <- cbind(valData,
                           obs =  dat[[colnamesResp]][testIDs],
                           pred = predict(modOut, datPreds[testIDs, ]))
          valData[, resid := pred - obs]

          ## Feature selection -- remove features (variables) with low SHAP values
          ## based on a quantile threshold
          shap_values <- shap.values(modOut, datPreds) |>
            Cache(omitArgs = formalArgs(shap.values),
                  .functionName = .functionNameHelper("shap.values", "xgboost", kFold),
                  .cacheExtra = c(dig, dig2, cols2keep))
          meanSHAP <- shap_values$mean_shap_score

          if (calcThresh) {
            SHAPthresh <- quantile(meanSHAP, prob = SHAPthresh)
            calcThresh <- FALSE ## we only calculate the threshold once
          }

          lowSHAPcols <- names(which(meanSHAP < SHAPthresh))
          if (length(lowSHAPcols)) {
            cols2keep <- setdiff(colnames(datPreds), lowSHAPcols)
            datPreds <- datPreds[, ..cols2keep]

            message("Removing features with SHAP < ", SHAPthresh, ":\n",
                    paste(lowSHAPcols, collapse = ", "))
          }

        }

        ## more outputs
        shapContrib <- shap_values$shap_score
        shapContrib <- shapContrib[, -"(Intercept)"]
        shap_long <- shap.prep(shap_contrib = shapContrib, X_train = datPreds) |>
          Cache(omitArgs = formalArgs(shap.prep),
                .functionName = .functionNameHelper("shap.prep", kFold),
                .cacheExtra = c(dig, dig2, cols2keep))

        list(valData = valData,
             mod = modOut,
             shap_values = shap_values,
             shap_long = shap_long)
      })
  )

  return(mm)
}

#' Calculate confidence matrices from continuous prediction from
#' XGBoost model
#'
#' @param mod a list containing `valData`, a `data.table` containing the
#'   validation data, predictions, observations and residuals from an xgboost
#'   cross-validation fold.
#' @param classMap a data.table of class to continuous value correspondences,
#'   with column names being `classVar` and `contVar`
#' @param classes vector of classes.
#' @param classVar the class variable/column name
#' @param contVAR the continuous variable/column name
#'
#' @returns
#' @export
#'
#' @examples
xgboostConfMat <- function(mod, classMap, classes, classVar = "SEV_CLASS", contVAR = "SEV_PROP") {
  predictionsDT <- copy(mod$valData)
  setnames(predictionsDT, "obs", "SEV_PROP")

  mappedClasses <- classMap[match(predictionsDT[[contVAR]], classMap[[contVAR]]), ..classVar]
  predictionsDT[, obsCLASS := mappedClasses]

  ## convert to classes, using the quantiles corresponding to the observed class proportions
  ## accumulate proportions to get probabilities
  quantProbs <- cumsum(table(predictionsDT$obsCLASS)/nrow(predictionsDT))
  classRanges <- c(0, quantile(predictionsDT$pred, probs = quantProbs))

  predictionsDT[, predCLASS := cut(pred, breaks = classRanges,
                                   include.lowest = TRUE, right = FALSE)]  ## classify as with intervals as ],]

  ## convert to numbered factor (subtracting one, because classes are 0-5)
  predictionsDT[, predCLASS := as.numeric(predCLASS)-1]
  predictionsDT[, `:=`(obsCLASS = factor(obsCLASS, levels = classes),
                       predCLASS = factor(predCLASS, levels = classes))]

  validMetrics <- caret::multiClassSummary(predictionsDT[, list(obs = obsCLASS,
                                                                pred = predCLASS)],
                                           lev = classes)
  ## calculate confusion matrix
  confMatrix <- caret::confusionMatrix(data = predictionsDT$predCLASS,
                                       reference = predictionsDT$obsCLASS)
}

#' Tune XGBoost parameters with `caret`
#'
#' Tuning is done in 3 steps.
#' Step 1. Tune learning rate (called `learning_rate` in `xgboost`
#' and `eta` in `caret`)
#' Step 2. Take the best learning rate and tune all other parameters
#' (from those passed to `xgboost` by `[caret::train()]`) except
#' `rnounds` and `sample` (which is always kept as 1).
#' Step 3. Take the best parameters values from Steps 1 and 2 and tune
#' `nrounds`.
#'
#' @returns a data.frame of best parameter values.
#'
#' @inheritParams runXGBOOST
#' @importFrom caret trainControl train caretTheme
#' @importFrom reproducible Cache
#' @importFrom lattice trellis.par.set
.tunexgboost <- function(dig, dat, colnamesResp, figDir) {
  ## use devtools::load_all("C:/Users/cbarros/GitHub/caret/pkg/caret/")
  ## bug reported at: https://github.com/topepo/caret/issues/1412

  savePlot <- FALSE
  if (!is.null(figDir)) {
    dir.create(figDir, showWarnings = FALSE, recursive = TRUE)
    savePlot <- TRUE
  }

  ## Step 1. tune learning rate.
  ## eta = learning rate.
  param_grid1 <- data.frame(nrounds = 200,
                            eta = seq(0.01, 1, by =  0.01),
                            ## defaults in xgboost:
                            max_depth = 6,
                            gamma = 0,
                            colsample_bytree = 1,
                            min_child_weight = 1,
                            subsample = 1)

  xgb_trcontrol <- trainControl(
    method = "cv",
    number = 5,
    verboseIter = TRUE,
    returnData = FALSE,
    returnResamp = "final",
    allowParallel = TRUE,
    savePredictions = "final"
  )
  message(cyan("Tuning learning rate..."))
  st <- system.time(
    {
      xgb_tuned <- train(SEV_PROP ~ .,
                         data = as.data.frame(dat),
                         trControl = xgb_trcontrol,
                         tuneGrid = param_grid1,
                         method = "xgbTree"
      ) |>
        Cache(omitArgs = c("data", "x"),
              .functionName = .functionNameHelper("train", "tune_learningrate"),
              .cacheExtra = c(dig),
              showSimilar = TRUE,
              cacheId = "8a518e3d96830586",   ## fix for now
              cacheSaveFormat = "rds")
    }
  )

  paramsF <- xgb_tuned$bestTune
  message(cyan("Finished in", st[["elapsed"]], "sec."))

  ## save tuning output
  if (savePlot) {
    png(file.path(figDir, "tuning_learningRate.png"), height = 4, width = 6,
        units = "in", res = 300)
    trellis.par.set(caretTheme())
    print(plot(xgb_tuned))
    dev.off()
  }

  ## Step 2. fix best learning rate and vary the rest
  browser()
  param_grid2 <- expand.grid(nrounds = 200,
                             max_depth = c(1:10),
                             eta = paramsF$eta,
                             gamma = c(0, 0.1, 1, 2),#, 5, 10), ## tested with more initially, but not necessary
                             colsample_bytree = c(0.1, 0.5, 1),
                             min_child_weight = c(0, 1, 2, 5),
                             subsample = 1)

  ## tune other parameters
  for (i in 1:3) gc(reset = TRUE)
  message(cyan("Tuning remaining XGBoost parameters..."))
  st <- system.time(
    {
      xgb_tuned <- train(SEV_PROP ~.,
                         data = as.data.frame(dat),
                         trControl = xgb_trcontrol,
                         tuneGrid = param_grid2,
                         method = "xgbTree"
      ) |>
        Cache(omitArgs = c("data", "x"),
              .functionName = .functionNameHelper("train", "tune_all"),
              .cacheExtra = c(dig),
              showSimilar = TRUE,
              ## cacheId = "d76ffa84709d8db0", # more params
              cacheId = "cd64c0826cfcecc2", #param grid above.
              cacheSaveFormat = "rds")
    }
  )

  paramsF <- xgb_tuned$bestTune
  message(cyan("Finished in", st[["elapsed"]], "sec."))   ## about 4hrs

  ## save tuning output
  if (savePlot) {
    png(file.path(figDir, "tuning_all.png"), height = 12, width = 12,
        units = "in", res = 300)
    trellis.par.set(caretTheme())
    print(plot(xgb_tuned))
    dev.off()
  }

  ## Step 3. vary only no. rounds
  param_grid3 <- expand.grid(nrounds = c(100, 200, 500, 1000, 1500),
                             max_depth = paramsF$max_depth,
                             eta = paramsF$eta,
                             gamma = paramsF$gamma,
                             colsample_bytree = paramsF$colsample_bytree,
                             min_child_weight = paramsF$min_child_weight,
                             subsample = 1)

  ## tune other parameters
  for (i in 1:3) gc(reset = TRUE)
  message(cyan("Tuning no. rounds (trees)..."))
  st <- system.time(
    {
      xgb_tuned <- train(SEV_PROP ~.,
                         data = as.data.frame(dat),
                         trControl = xgb_trcontrol,
                         tuneGrid = param_grid3,
                         method = "xgbTree"
      ) |>
        Cache(omitArgs = c("data", "x"),
              .functionName = .functionNameHelper("train", "tune_nrounds"),
              .cacheExtra = c(dig),
              showSimilar = TRUE,
              cacheId = "3499577823e6a19d",
              cacheSaveFormat = "rds")
    }
  )

  paramsF <- xgb_tuned$bestTune
  message(cyan("Finished in", st[["elapsed"]], "sec."))
  message(cyan("Best parameters:"))
  message(cyan(paste0(capture.output(paramsF), collapse = "\n")))

  ## save tuning output
  if (savePlot) {
    png(file.path(figDir, "tuning_nrounds.png"), height = 4, width = 6,
        units = "in", res = 300)
    trellis.par.set(caretTheme())
    print(plot(xgb_tuned))
    dev.off()
  }

  for (i in 1:3) gc(reset = TRUE)
  return(paramsF)
}

## GPBoost wrapper -----------------------
#' Wrapper for gpboost
#' @inheritParams runXGBOOST
runGPBOOST <- function(dat, dig, nFolds = 5, colnamesPred, colnamesResp = "SURV_PROP",
                       colnamesGrp, colnamesCat, SHAPthresh = 0) {

  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  # Setup k-folds
  indexNames <- c("allData", "evalData")

  ## create folds and make a list with indices of full dataset and eahc fold
  trainIndexK <- createFolds(dat[[colnamesResp]], k = nFolds, list = TRUE, returnTrain = FALSE)
  trainIndexK <- Map(tr = trainIndexK, function(tr) {
    list(seq(NROW(dat)), tr) |> setNames(indexNames)
  })

  ## sample columns after setting seed for caching (if different, then cache is triggered)
  colOrder <- colnames(dat)
  colOrder <- sample(colOrder)
  dat <- dat[, ..colOrder]
  dig <- .robustDigest(dat)

  ## subset predictor data
  if (missing(colnamesPred)) {
    message(cyan("No predictor columns ('colnamesPred') provided.",
                 "\nAssuming all but 'colnamesResp' and 'colnamesGrp' are predictors"))
    colnamesPred <- setdiff(c(colnames(dat), colnamesGrp), colnamesResp)
  }
  datPreds <- dat[, ..colnamesPred]

  st <- system.time(
    mm <- pmap(
      list(dataFolds = trainIndexK, kFold = seq(nFolds)),
      function(dataFolds, kFold) {
        ## get row IDs for training data (allData) and testing data
        allDataIDs <- dataFolds[[indexNames[[1]]]]
        testIDs <- dataFolds[[indexNames[[2]]]]   ## eval data
        dig2 <- .robustDigest(dataFolds)

        lowSHAPcols <- 1
        calcThresh <- TRUE
        # SHAPthresh <- 0.25  ## test
        cols2keep <- colnames(datPreds)
        modelGPBfit <- NULL

        while (length(lowSHAPcols)) {
          ## TODO: export RMSEs; go back to previous model if AUC decreases after removing features
          browser()

          gp_train <- gpb.Dataset(as.matrix(datPreds[allDataIDs]),
                                  label = dat[[colnamesResp]][allDataIDs],
                                  categorical_feature = colnamesCat,
                                  free_raw_data = FALSE)
          gp_model <- GPModel(group_data = dat[[colnamesGrp]][allDataIDs],
                              likelihood = "zero_inflated_gamma",
                              free_raw_data = FALSE)
          gp_model$set_optim_params(params = list(trace = TRUE))

          ## TODO: do a 2 step approach. tun learning_rate first, then the rest.
          ## define a parameter grid and fixed parameters to optimize

          N <- length(dat[[colnamesResp]][allDataIDs])
          param_grid <- list("learning_rate" = c(0.01, 0.1, 1),
                             "min_data_in_leaf" = c(1, 10, 100, 1000),
                             "max_depth" = c(-1, 1:5),
                             "num_leaves" = 2^(1:10),
                             "lambda_l2" = c(0, 1, 10, 100),
                             "max_bin" = c(250, 500, 1000, min(N,10000)),
                             "line_search_step_length" = c(TRUE, FALSE))

          ## TODO: RE-RUN AFTER PRESENTATION
          opt_params <- loadFromCache(cacheId = "da3593ff30bbc3f7")
          # opt_params <- gpb.grid.search.tune.parameters(param_grid = param_grid,
          #                                               data = gp_train,
          #                                               gp_model = gp_model,
          #                                               num_try_random = 1,
          #                                               folds = list(testIDs),
          #                                               nrounds = 20,
          #                                               early_stopping_rounds = 10,
          #                                               verbose_eval = 1,
          #                                               metric = "rmse") |>
          #   Cache(omitArgs = c("data", "gp_model"),
          #         .functionName = .functionNameHelper("gpb.tune.parameters"),
          #         .cacheExtra = c(dig, dig2, cols2keep, colnamesGrp, kFold),
          #         showSimilar = TRUE,
          #         cacheSaveFormat = "rds")

          ## Stage 1: run cross-validation to
          ## (i) determine to optimal number of iterations
          params <- opt_params$best_params
          modelGPBcv <- gpb.cv2(data = gp_train
                                , gp_model = gp_model
                                , metric = c("auc", "rmse")
                                , eval = c("auc", "rmse")
                                , nrounds = 200
                                , folds = list(testIDs)
                                , early_stopping_rounds = opt_params$best_iter + 10
                                , params = params
                                # , use_gp_model_for_validation = FALSE ## this causes a fatal error
          ) |>
            Cache(omitArgs = c("data", "gp_model"),
                  .functionName = .functionNameHelper("gpb.cv2", kFold),
                  .cacheExtra = c(dig, dig2, cols2keep, colnamesGrp),
                  showSimilar = TRUE,
                  cacheSaveFormat = "rds")

          print(paste0("Optimal number of iterations: ", modelGPBcv$best_iter,
                       ", best test error: ", modelGPBcv$best_score))

          ## Stage 2: Train tree-boosting model
          modelGPBfit2 <- (function(){
            out <- gpb.train(data = gp_train,
                             gp_model = gp_model,
                             nrounds = modelGPBcv$best_iter + 10,
                             params = params)

            list(modelGPBfit = out, gp_model = gp_model)
          })() |>
            Cache(.functionName = .functionNameHelper("gpb.train", kFold),
                  .cacheExtra = c(dig, dig2, cols2keep),
                  showSimilar = TRUE,
                  cacheSaveFormat = "rds")


          if (is.null(modelGPBfit)) {
            modelGPBfit <- modelGPBfit2
            AUCout <- tail(attr(modelGPBfit, "evaluation_log"), 1)$train_auc
          }

          ## get last AUC
          AUCout2 <- tail(attr(modelGPBfit2, "evaluation_log"), 1)$train_auc
          if (AUCout2 < AUCout) {
            message("AUC decreased after removing features.\n",
                    "  The previous model will be retained, instead")
          } else {
            modelGPBfit <- modelGPBfit2
          }
          AUCout <- AUCout2

          ## calculate predictions and residuals
          valData <- datPreds[testIDs,]
          ## here: error in predict
          pred <- modelGPBfit$modelGPBfit$predict(
            data = as.matrix(datPreds[testIDs, ]),
            group_data_pred = dat[[colnamesGrp]][testIDs]
          )
          valData <- cbind(valData,
                           obs =  dat[[colnamesResp]][testIDs],
                           pred = pred$response_mean
          )
          valData[, resid := pred - obs]

          ## Feature selection -- remove features (variables) with low SHAP values
          ## based on a quantile threshold
          shap_values <- shap.values(modelGPBfit$modelGPBfit, as.matrix(datPreds[allDataIDs])) |>
            Cache(omitArgs = formalArgs(shap.values),
                  .functionName = .functionNameHelper("shap.values", "gpboost", kFold),
                  .cacheExtra = c(dig, dig2, cols2keep))
          meanSHAP <- shap_values$mean_shap_score

          if (calcThresh) {
            SHAPthresh <- quantile(meanSHAP, prob = SHAPthresh)
            calcThresh <- FALSE ## we only calculate the threshold once
          }

          lowSHAPcols <- names(which(meanSHAP < SHAPthresh))
          if (length(lowSHAPcols)) {
            cols2keep <- setdiff(colnames(datPreds), lowSHAPcols)
            datPreds <- datPreds[, ..cols2keep]
            colnamesCat <- intersect(colnamesCat, colnames(datPreds))

            message("Removing features with SHAP < ", SHAPthresh, ":\n",
                    paste(lowSHAPcols, collapse = ", "))
          }
        }

        ## more outputs
        shapContrib <- shap_values$shap_score
        shapContrib <- suppressWarnings(shapContrib[, -"(Intercept)"])
        shap_long <- shap.prep(shap_contrib = shapContrib, X_train = datPreds[allDataIDs]) |>
          Cache(omitArgs = formalArgs(shap.prep),
                .functionName = .functionNameHelper("shap.prep", "gpboost", kFold),
                .cacheExtra = c(dig, dig2, cols2keep))

        list(valData = valData,
             modelGPBfit = modelGPBfit,
             modelGPBcv = modelGPBcv,
             shap_values = shap_values,
             shap_long = shap_long)
      })
  )

  return(mm)
}

.predfunGAMLSS <- function(x, newdata) {
  ## note shapr::explain documentation is wrong, first arg needs to be called x
  dig <- .robustDigest(newdata)
  dig2 <- .robustDigest(x)
  preds <- predictAll(x, newdata = newdata, type = "response") |>
    Cache(omitArgs = formalArgs(predictAll),
          .cacheExtra = c(dig, dig2))

  calcMeanBEINF(preds$mu, preds$nu, preds$tau)
}

.modelspecsfunGAMLSS <- function(x) {
  featlabs <- c(labels(x$mu.terms),
                labels(x$nu.terms),
                labels(x$sigma.terms),
                labels(x$tau.terms))
  featlabs <- unlist(strsplit(featlabs, split = ":"))
  featlabs <- unique(featlabs)

  featclass <- c(attr(x$mu.terms, "dataClasses")[-1],
                 attr(x$nu.terms, "dataClasses")[-1],
                 attr(x$sigma.terms, "dataClasses")[-1],
                 attr(x$tau.terms, "dataClasses")[-1])
  featclass <- featclass[featlabs]

  featlabs <- sub(".*(FIRE_NAME).*", "\\1",  featlabs)
  names(featclass) <- sub(".*(FIRE_NAME).*", "\\1",  names(featclass))
  featclass["FIRE_NAME"] <- "factor"

  feature_specs <- list()
  feature_specs$labels <- featlabs
  feature_specs$classes <- featclass


  m <- length(feature_specs$labels)
  feature_specs$factor_levels <- setNames(vector("list", m), feature_specs$labels)
  feature_specs$factor_levels[feature_specs$classes == "factor"] <- NA # model object doesn't contain factor levels info

  return(feature_specs)
}

.functionNameHelper <- function(..., sep = "_") {
  paste(..., sep = sep)
}

#' Calculate confidence matrices from continuous prediction from
#' GPBoost model
#'
#' @param mod a list containing `valData`, a `data.table` containing the
#'   validation data, predictions, observations and residuals from an xgboost
#'   cross-validation fold.
#' @param classMap a data.table of class to continuous value correspondences,
#'   with column names being `classVar` and `contVar`
#' @param classes vector of classes.
#' @param classVar the class variable/column name
#' @param contVAR the continuous variable/column name
#'
#' @returns
#' @export
#'
#' @examples
gpboostConfMat <- function(mod, classMap, classes, classVar = "SEV_CLASS", contVAR = "SEV_PROP") {
  predictionsDT <- copy(mod$valData)
  predictionsDT[, `:=`(pred = 1 - pred,
                       obs = 1 - obs)]  ## back transform to severity
  setnames(predictionsDT, "obs", contVAR)

  mappedClasses <- classMap[match(predictionsDT[[contVAR]], classMap[[contVAR]]), ..classVar]
  predictionsDT[, obsCLASS := mappedClasses]

  ## convert to classes, using the quantiles corresponding to the observed class proportions
  ## accumulate proportions to get probabilities
  quantProbs <- cumsum(table(predictionsDT$obsCLASS)/nrow(predictionsDT))
  classRanges <- c(0, quantile(predictionsDT$pred, probs = quantProbs))

  predictionsDT[, predCLASS := cut(pred, breaks = classRanges,
                                   include.lowest = TRUE, right = FALSE)]  ## classify as with intervals as ],]

  ## convert to numbered factor (subtracting one, because classes are 0-5)
  predictionsDT[, predCLASS := as.numeric(predCLASS)-1]
  predictionsDT[, `:=`(obsCLASS = factor(obsCLASS, levels = classes),
                       predCLASS = factor(predCLASS, levels = classes))]

  validMetrics <- caret::multiClassSummary(predictionsDT[, list(obs = obsCLASS,
                                                                pred = predCLASS)],
                                           lev = classes)
  ## calculate confusion matrix
  confMatrix <- caret::confusionMatrix(data = predictionsDT$predCLASS,
                                       reference = predictionsDT$obsCLASS)

  list(validMetrics = validMetrics, confMatrix = confMatrix)
}



## Plotting functions -------------
plotFun <- function(firePoints, var, fireID, varTitle) {
  fireID <- sub("_(obs|pred|p1|varY)", "", fireID)
  ggplot(firePoints[firePoints$FIRE_NAME == fireID]) +
    geom_spatvector(aes_string(colour = var)) +
    scale_color_distiller(palette = "RdYlBu", direction = -1,
                          limits = c(0,1)) +
    theme(axis.text.x = element_text(hjust = 1, angle = 20)) +
    labs(colour = "", x = "Longitude", y = "Latitude",
         title = varTitle)
}


plotResiduals <- function(valData, filename = NULL) {
  plot1 <- ggplot(valData, aes(x = pred, y = resid)) +
    geom_point() +
    theme_bw()

  plot2 <- ggplot(valData, aes(x = 1:nrow(valData), y = resid)) +
    geom_point() +
    theme_bw()

  plot3 <- ggplot(valData, aes(x = resid)) +
    geom_density() +
    theme_bw()

  plot4 <- ggplot(valData, aes(sample = resid)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw()

  if (!is.null(filename)) {
    png(filename, width = 7, height = 7, units = "in", res = 300)
    on.exit(dev.off(), add = TRUE)
  }
  ggarrange(plotlist = list(plot1, plot2, plot3, plot4))
}
