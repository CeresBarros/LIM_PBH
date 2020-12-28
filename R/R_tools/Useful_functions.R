## ------------------------------------------------------
## USEFUL FUNCTIONS
##
## Ceres: March 2020
## ------------------------------------------------------

## this script should be sourced
require(raster); require(tools)
require(sp); require(dismo)
require(sf); require(dplyr)

## STANDARDIZE DATA TO 0-1 RANGE ------------------------
## x is a numeric vector
## NOTE: this function is available on amc.
rescale01 <- function(x) {
  xx <- (x - min(x))/(max(x) - min(x))
  return(xx)
}

## CHECK PROJECTIONS ------------------------------------
## sfObj.list is a list of spatial objects
checkProjections <- function(sfObj.list){
  projs <- sapply(sfObj.list, FUN = function(x) {
    return(
      eval(expr = parse(text = paste0("projection(", x, ")")))
    )
  })
  return(projs)
}


## CROP & MASK TO STUDY AREA ----------------------------
## study.area is a "Raster*" or "Spatial*" object
## tocrop must be a Raster
## method is passed to projectRaster - might need to be changed for factors
cropToStudyArea <- function(study.area, tocrop, method = "bilinear") {
  temp <- tocrop

  if(class(study.area) == class(temp) & class(study.area) == "RasterLayer") {
    temp <- projectRaster(from = temp, to = study.area, crs = crs(study.area), method = method)
  } else {
    if(crs(study.area)@projargs != crs(tocrop)@projargs){
      temp <- projectRaster(tocrop, crs = crs(study.area))
    }
  }
  temp <- crop(x = temp, y = study.area)
  temp <- mask(x = temp, mask = study.area)

  return(temp)
}

## RASTER TO BINARY MATRIX ------------------------------
## converts a vector of values into a binary "presence/absence" matrix
## x is a vector
vector2binmatrix <- function(x) {
  x <- as.character(x)
  return(model.matrix( ~ x-1))
}


## FIND NEIGHBOURS IN MATRIX ----------------------------
## finds the 8 neighbours of each cell in a matrix
## output is a matrix of 8 rows, with a columns per cell of the input matrix,
## which is treated by columns

neighboursMatrix = function(mat) {
  mat2 <- cbind(NA, rbind(NA, mat ,NA), NA)   ## makes a border of NAs
  addresses <- expand.grid(x = 1:nrow(mat), y = 1:ncol(mat)) ## all matrix coordinates
  neighs <- c()
  for(i in 1:-1) {
    for(j in 1:-1) {
      if(i != 0 || j != 0) {
        neighs <- rbind(neighs,mat2[addresses$x+i+1+nrow(mat2)*(addresses$y+j)])   ## each column contains the neighboors a cell (going by columns in mat)
      }
    }
  }
  return(neighs)
}



## JOIN SPATIAL OBJECTS FUNCTION -----------------------
## joins shapefiles or rasters
## files is a character string of file names (does not require extension)
## folder is the folder where they can be found

loadBindSpatialObjs <- function(files, destinationPath, urls = NULL) {
  ## name URLs with file names
  if(length(urls) > 1) {
    names(urls) = files

    if(all(grepl(".shp", files))) {
      sfObj.ls <- lapply(files, FUN = function(targetFile) {
        prepInputs(targetFile = file.path(destinationPath, targetFile),
                   url = urls[targetFile], destinationPath = destinationPath,
                   fun = "shapefile", pkg = "raster")
      })
      if(length(unique(sapply(sfObj.ls, FUN = function(sfObj) as.character(crs(sfObj))))) == 1) {
        joined = do.call(bind, sfObj.ls)
      } else(stop("Files do not share the same projection"))

    }  else {
      ## check if all are raster files
      if(all(grepl(".grd", files)) |
         all(grepl(".asc", files)) |
         all(grepl(".tif", files)) |
         all(grepl(".img", files))) {

        sfObj.ls <- lapply(files, FUN = function(targetFile) {
          prepInputs(targetFile = file.path(destinationPath, targetFile),
                     url = urls[targetFile], destinationPath = destinationPath,
                     fun = "raster", pkg = "raster")
        })

        ## check projections match and do the join
        if(length(unique(sapply(sfObj.ls, FUN = function(sfObj) as.character(crs(sfObj))))) == 1) {
          joined = do.call(bind, sfObj.ls)
        } else(stop("Files do not share the same projection"))

      } else stop("All files should be in the same format (.shp, .grd, .asc, .tif or .img)")
    }

  } else {
    if(all(grepl(".shp", files))) {
      sfObj.ls <- lapply(files, FUN = function(targetFile) {
        prepInputs(targetFile = targetFile,
                   url = urls, destinationPath = destinationPath,
                   fun = "shapefile", pkg = "raster")
      })
      if(length(unique(sapply(sfObj.ls, FUN = function(sfObj) as.character(crs(sfObj))))) == 1) {
        joined = do.call(bind, sfObj.ls)
      } else(stop("Files do not share the same projection"))

    }  else {
      ## check if all are raster files
      if(all(grepl(".grd", files)) |
         all(grepl(".asc", files)) |
         all(grepl(".tif", files)) |
         all(grepl(".img", files))) {

        sfObj.ls <- lapply(files, FUN = function(targetFile) {
          prepInputs(targetFile = targetFile,
                     url = urls, destinationPath = destinationPath,
                     fun = "raster", pkg = "raster")
        })

        ## check projections match and do the join
        if(length(unique(sapply(sfObj.ls, FUN = function(sfObj) as.character(crs(sfObj))))) == 1) {
          joined = do.call(bind, sfObj.ls)
        } else(stop("Files do not share the same projection"))

      } else stop("All files should be in the same format (.shp, .grd, .asc, .tif or .img)")
    }
  }
  return(joined)
}


## DRAW CONVEX HILL AROUND POLYGONS -----------------------
## draws a convex hull aroung vertice points of a polygon shape file.
## x must be a SpatialPolygons, or SpatialPolygonsDF
## NOTE: this function as been passed to amc.
outerBuffer <- function(x) {
  if(class(x) == "SpatialPolygons" | class(x) == "SpatialPolygonsDataFrame") {
    ## Get polygon vertices
    pts <- SpatialPoints(do.call(rbind, lapply(x@polygons, FUN = function(x) {
      return(x@Polygons[[1]]@coords)
    })))

    ## Draw convex hull around points and extract polygons slot
    hull <- polygons(convHull(pts))

    return(hull)
  } else(stop("x must be a SpatialPolygons, or SpatialPolygonsDF"))
}


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


## CLEAN SF OBJECT FROM DUPLICATED COLUMNS ------------------------
## this fucnton removes potentially duplicated data columns in an sf object
## sfObj is an sf object from the sf package
sfRmDupNACols <- function(sfObj) {
  transposed <- t(st_set_geometry(sfObj, NULL))
  ## get duplicates and NAs (note that cols became rows)
  dupCols <- duplicated(transposed)
  NAcols <- sapply(sfObj[,, drop = TRUE], FUN = function(var) all(is.na(var)))

  if (any(dupCols)) {
    transposed <- transposed[!dupCols,]
    dataSF <- data.frame(t(transposed), stringsAsFactors = FALSE)

    ## convert columns to numeric where appropriate
    numCols <- intersect(names(dataSF), names(which(sapply(st_set_geometry(sfObj, NULL), is.numeric))))
    for(col in numCols) {
      dataSF[, col] <- as.numeric(dataSF[, col])
    }

    ## check for NA cols that were not removed
    NAcols <- sapply(dataSF, FUN = function(var) all(is.na(var)))
    if (any(NAcols))
      dataSF <- dataSF[, !NAcols]

    ## add geometry to data.frame to re-make sf object
    st_geometry(dataSF) <- st_geometry(sfObj)
    sfObj <- dataSF
  } else {
    if (any(NAcols)) {
      sfObj <- sfObj[, !NAcols]
    } else message("no duplicated or NA columns were found.")
  }
  return(sfObj)
}


## RENAME AND SUBSET FIELDS IN SF OBJECT ACCORDING TO TABLE------------------------
## sfObj is an sf object from the sf package
## namesTable is a two-column data.frame with the names to be replaced (1st column) and the new names (2nd column)
##        if there are NAs in the second columns, the original columns will be removed if rmNAs = TRUE
## rmNAs determines whether columns wiith missing new names are removed or not.

renameCleanSfFields <- function(sfObj, namesTable, rmNAs = TRUE) {
  origNames <- names(st_set_geometry(sfObj, NULL))
  if (!all(origNames %in% namesTable$Name_shp))
    stop("Some names in 'sfObj' are missing from 'namesTable'")

  ## get new names
  rownames(namesTable) <- namesTable[, 1]
  newNames <- as.character(namesTable[origNames, 2])

  ## if there are missing new names, either remove them, or maintain original name
  if (any(is.na(newNames))) {
    if (rmNAs) {
      toRM <-  origNames[is.na(newNames)]
      sfObj[, toRM] <- NULL

      ## re-do steps above
      origNames <- names(st_set_geometry(sfObj, NULL))
      newNames <- as.character(namesTable[origNames, 2])
    } else {
      newNames[is.na(newNames)] <- origNames[is.na(newNames)]
    }
  }

  names(sfObj)[names(sfObj) %in% origNames] <- newNames

  return(sfObj)
}


## VALIDATE GEOMETRIES ------------------------
## sfObj is an sf object from the sf package
## dim is used for faster caching
validateGeomsSf <- function(sfObj, dim) {
  ## Any corrupt or invalid geometries?
  if (any(is.na(st_is_valid(sfObj))) |
      any(na.omit(st_is_valid(sfObj)) == FALSE)) {
    message("Invalid gemoetries found. Attempting to make valid using st_make_valid")
    sfObj <- st_make_valid(sfObj)
  }

  if (any(is.na(st_is_valid(sfObj))) |
      any(na.omit(st_is_valid(sfObj)) == FALSE))
    stop("st_make_valid did not work, please check what's wrong.")

  sfObj
}

## DOWNLOAD KMZ AND CONVERT TO SHAPEFILE --------------------
## url is a google drive URl
## archive is the .zip file name where the .kmz is contained
## destination path is the path to where archive will downloaded to, and where the final shapefile will be saved.
prepKMZ2shapefile <- function(url, archive, destinationPath) {
  ## check
  if (class(url) != "character" | is.null(url))
    stop("Provide url as a character string")
  if (!grepl("\\.zip$", archive))
    stop("archive should be the name of the .zip file in url")
  if (is.null(destinationPath)) {
    message("archive will be downloaded to tempdir()")
    destinationPath <- tempdir()
  }

  archivePath <- file.path(destinationPath, archive)
  downloaded_file <- googledrive::drive_download(file = googledrive::as_id(url),
                                                 path = archivePath,
                                                 overwrite = TRUE)
  if (downloaded_file$local_path != archivePath |
      downloaded_file$name != archive)
    stop(paste0("Downloaded file name (", downloaded_file$name, ") and archive (",
                archive, ") do not match."))

  ## unzip .zip and .kmz
  fileKMZ <- unzip(archivePath, exdir = destinationPath)
  fileKML <- unzip(fileKMZ, exdir = destinationPath)

  ## load as sf
  sfObj <- st_read(fileKML)
  if (any(names(sfObj) %in% "Description"))
    sfObj$Description <- NULL  ## weird unnecessary column
  ## convert to shapefile
  shpObj <- as(st_zm(sfObj), "Spatial")

  ## delete unecessary .kmz/.kml files and save .shp
  file.remove(fileKML, fileKMZ)
  shpFile <- sub(".zip", ".shp", archive, fixed = TRUE)
  shapefile(shpObj, filename = file.path(destinationPath, shpFile), overwrite = TRUE)

  return(shpObj)
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
      plan(multiprocess(gc = TRUE))
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
      plan(multiprocess(gc = TRUE), workers = cores)
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
    ngbSEVDT <- rbindlist(ngbSEVList, use.names = TRUE, fill = TRUE)
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
                               fire = fireID)
  ngbhoodBurnsDT <- na.omit(ngbhoodBurnsDT)
  setnames(ngbhoodBurnsDT, old = c("ngbPropBurns", "fire"),
           new = c(paste0("ngbPropBurns", sevColID, "_", dist, "m"),
                   fireColID))

  ## checks
  if (length(unique(ngbhoodBurnsDT$pixID)) != length(i)) {
    noMissing <- sum(!sevPoints$pixID[i] %in% ngbhoodBurnsDT$pixID)
    warning(paste(noMissing, "points in fire perimeter were not converted to raster pixels.",
                  "\nThis is probably due to more than one point falling in the same cell"))
  }
  message(paste0("Done!"))
  return(ngbhoodBurnsDT)
}
