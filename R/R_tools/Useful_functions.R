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




# importFrom caret createFolds
# importFrom purrr pmap
# importFrom pROC roc
runXGBOOST <- function(dat, dig, nFolds = 5, colnamesPred, colnamesResp = "SEV_PROP"
                       #, crossValType = c("time-ordered", "crossValidation")  ## automatically decided
) {
  # dat <- dat[sample(NROW(dat), size = 1e4)]
  # tt <- table(dat$SEV_PROP)

  # Add dummy variables for factor columns -- i.e., the random effects
  if (all(sapply(dat, is.numeric)) %in% FALSE)
    dat <- #{
      model.matrix(~ . + 0, data = dat) |>
      Cache(omitArgs = c("object", "data", "x"), .cacheExtra = dig) # Creates dummy variables

  dat <- as.data.table(dat)

  savedSeed <- .Random.seed
  on.exit(assign(".Random.seed", savedSeed, envir = .GlobalEnv), add = TRUE)
  set.seed(12345) # so kfolds are same, so Caching works correctly below; if dat3 changes number of rows,
  # it will be a totally different sequence; but it will be the same sequence
  # if number of rows doesn't change

  # Setup k-folds
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
  colOrder <- setdiff(colnames(dat), c("pixID", yearColname))
  colOrder <- sample(colOrder)
  dat <- dat[, ..colOrder]
  dig <- .robustDigest(dat)

  ## subset predictor data
  if (missing(colnamesPred)) {
    message(cyan("No predictor columns ('colnamesPred') provided.",
                 "\nAssuming all but 'colnamesResp' are predictors"))
    colnamesPred <- setdiff(colnames(dat), colnamesResp)
  }
  datPreds <- dat[, ..colnamesPred]

  st <- system.time(
    mm <- pmap(
      list(dataFolds = trainIndexK, kFold = seq(nFolds)),
      function(dataFolds, kFold)#, indexName = indexNames, digInner = dig, dat3ForxgboostInner = dat,
        #         datPredsInner = datPreds)
      {
        ## get row IDs for training data (allData) and testing data
        allDataIDs <- dataFolds[[indexNames[[1]]]]
        testIDs <- dataFolds[[indexNames[[2]]]]   ## eval data
        digValInd <- .robustDigest(dataFolds) ## TODO ask Eliot: should be eval set?

        # xgboost objects do not save with `qs` ... must be `rds`
        # opt <- options(reproducible.cacheSaveFormat = "rds")
        # on.exit(options(opt)) # redundant; but necessary if it fails during fit
        mMSE <- xgboost(x = datPreds[allDataIDs],
                        , y = dat[, get(colnamesResp)][allDataIDs]
                        # objective = "count:poisson", #
                        # objective = "reg:tweedie",
                        , nthread = 10
                        , eval_set = testIDs, # 0.2, # should be evalData
                        , monitor_training = TRUE
                        # verbosity = 0,
                        , eval_metric = c("auc", "logloss")
                        , nrounds = 100
                        , max_depth = 2
                        # reg_lambda = 0.5,
                        # learning_rate = 0.15
        ) #|> Cache(omitArgs = c("x", "y", "eval_set"),
                   # # ) |> list()} |> Cache(omitArgs = c("x", "y", "eval_set", "..."),
                   # .functionName = functionNameHelper("xgboost", type, kFold),
                   # .cacheExtra = c(dig, digValInd, type),
                   # cacheSaveFormat = "rds")
                   #

        ## calculate predictions and residuals
        valData <- dat[testIDs,]
        valData <- cbind(valData,
                         obs =  dat[, get(colnamesResp)][testIDs],
                         pred = predict(mMSE, dat[testIDs, ]))
        valData[, resid := pred - get(colnamesResp)]

        ## calculate R^2 -- probably not adequate
        # sstot <- sum((valData$pred - mean(valData$obs))^2)
        # ssresid <- sum(valData$resid^2)
        # sprintf("percent variance explained, R^2: %1.1f%%", 100 * (1 - ssresid / sstot))

        # mMSE <- mMSE[[1]] # remove the list
        # Predict probabilities

        if (FALSE) {
          ignZeroAndOnes <- pmin(1, valData[, ignitions])
          # (roc_curvePoisson <- roc(ignZeroAndOnes, d$valData$predPoisson))
          (roc_curveTweedie <- roc(ignZeroAndOnes, valData[["predTweedie"]]))
        }
browser()
        shap_values <- shap.values(mMSE, datPreds) #|>
          # Cache(omitArgs = formalArgs(shap.values),
                # .functionName = functionNameHelper("shap.values", type, kFold),
                # .cacheExtra = c(dig, digValInd, type))
        shapContrib <- shap_values$shap_score
        shapContrib <- shapContrib[, -"(Intercept)"]
        shap_long <- shap.prep(# xgb_model = mMSE,
          shap_contrib = shapContrib, X_train = datPreds)  #|>
          # Cache(omitArgs = formalArgs(shap.prep),
                # .functionName = functionNameHelper("shap.prep", type, kFold),
                # .cacheExtra = c(dig, digValInd, type))
        # list(valData = valData, mod = mMSE, shap_long = shap_long)

        shap.plot.summary(shap_long, dilute = 20, scientific = TRUE)
        shap.plot.dependence(data_long = shap_long,
                             x = 'CMDsm',
                             y = 'FlamConif',
                             color_feature = 'Column_WV') +
                ggtitle("SHAP values of Pice_var vs. CMDsm")
      })
  )

  if (FALSE)  {
    # ggIgnit <-
    ggplot(df) +
      aes(x = value, y = predictedProb, group = variable, col = variable) +
      # ggIgnit <- ggplot(df) + aes(x = value, y = predictedProb, group = variable, col = variable) +
      geom_smooth(span = 1) +
      # ggplot2::geom_jitter(height = 0.00002, size = 0.5) +
      geom_point(size = 0.3) +
      # geom_line() +
      theme_bw() +
      scale_color_brewer(palette="Paired") +
      xlab(paste0("Scaled, centred")) +
      ylab(paste0("Predicted probability of ", type)) +
      guides(col = guide_legend(reverse = FALSE, col = factor(df$predictedProb)))
  }

  rocs <- lapply(mm, function(d) {
    ignZeroAndOnes <- pmin(1L, d$valData[[ignOrEscapeColName]])
    # (roc_curvePoisson <- roc(ignZeroAndOnes, d$valData$predPoisson))
    (roc_curveTweedie <- roc(ignZeroAndOnes, d$valData[["predTweedie"]]))
    roc_curveTweedie
  })

  tweedie <- mapply(r = rocs, function(r) {
    as.numeric(r$auc)
  })
  # poiss <- mapply(r = rocs, function(r) {
  #   as.numeric(r$roc_curvePoisson$auc)
  # })
  print(paste("mean roc: ", format(mean(tweedie), digits = 3)))
  # mean(poiss)
  mm2 <- Map(m = mm, function(m) m$mod)
  mm2 <- append(mm2, list(rocs = rocs))

  return(mm2)
}
