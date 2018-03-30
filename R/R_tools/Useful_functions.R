## ------------------------------------------------------
## USEFUL FUNCTIONS
##
## Ceres: Dec 2017
## ------------------------------------------------------

## this script should be sourced

## CHECK PROJECTIONS ------------------------------------
## obj.list is a list of spatial objects
checkProjections <- function(obj.list){
  projs <- sapply(obj.list, FUN = function(x) {
    return(
      eval(expr = parse(text = paste0("projection(", x, ")")))
    )
  })
  return(projs)
}


## CROP & MASK TO STUDY AREA ----------------------------
## study.area and tocrop are "Raster*" or "Spatial*" objects
cropToStudyArea <- function(study.area, tocrop) {
  require(raster)
  temp <- crop(x = tocrop, y = study.area)
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

loadBindSpatialObjs <- function(files, folder) {
  require(raster); require(tools)
  if(all(files %in% file_path_sans_ext(list.files(folder)))) {
    files2 <- grep(paste(files, collapse = "|"), list.files(folder), value = TRUE)
    
    ## check if all are shapefiles
    if(all(paste0(files, ".shp") %in% files2)) {
      files2 <- files2[files2  %in% paste0(files, ".shp")]
      
      ## load files
      obj.ls <- lapply(files2, FUN = function(f){
        eval(parse(text = paste0(
          file_path_sans_ext(f), " <- shapefile('", file.path(folder, f),"')"
        )))
        
        eval(parse(text = paste0("return(", file_path_sans_ext(f), ")")))
      })
      
      ## check projections match and do the join
      if(length(unique(sapply(obj.ls, FUN = function(obj) as.character(crs(obj))))) == 1) {
        joined = do.call(bind, obj.ls)
      } else(stop("Files do not share the same projection"))
      
    } else {
      ## check if all are raster files
      if(all(paste0(files, ".grd") %in% files2) |
         all(paste0(files, ".asc") %in% files2) |
         all(paste0(files, ".tif") %in% files2) |
         all(paste0(files, ".img") %in% files2)) {
        
        obj.ls <- lapply(files2, FUN = function(f){
          eval(parse(text = paste0(
            file_path_sans_ext(f), " <- raster('", file.path(folder, f),"')"
          )))
          
          eval(parse(text = paste0("return(", file_path_sans_ext(f), ")")))
        })
        
        ## check projections match and do the join
        if(length(unique(sapply(obj.ls, FUN = function(obj) as.character(crs(obj))))) == 1) {
          joined = do.call(bind, obj.ls)
        } else(stop("Files do not share the same projection"))
        
      } else stop("All files should be in the same format (.shp, .grd, .asc, .tif or .img)")
    }
    
    return(joined)
  } else(stop(paste("Can't find", 
                    files[files %in% file_path_sans_ext(list.files(folder))], 
                    "in", folder)))
}


## DRAW CONVEX HILL AROUND POLYGONS -----------------------
## draws a convex hull aroung vertice points of a polygon shape file.
## x must be a SpatialPolygons, or SpatialPolygonsDF
outerBuffer <- function(x) {
  require(sp); require(dismo)
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
## Method from Andison (2012), defines a fire event composed of
## disturbed patches (severity/mortality >= 95%)
## island remnants (severity < 95%, surrounded by disturbed patches)
## matrix remnants (undisturbed, and partially surrouned by disturbed patches)

## sf.obj is a Simple Features object from sf package
## fireNAMES should be a column/attribute with fire ID/names
## crsProj is a string of the projection to use. If NULL (default) will use the same projection as sf.obj, otherwise sf.obj will be reprojected
## buff.dist if the buffer distance to define fire events
## PLOT, SAVE and overwrite determine if plotting, saving and overwriting should be done
## outputDIR and fileNAME define the directory and fileNAME to save the fire events shapefile (".shp" will be added to the name string),
## outputDIR will be created if non-existent

defineFireEvents <- function(sf.obj, fireNAMES = NULL, fireVARS = NULL,crsProj = NULL, buff.dist = NULL, PLOT = TRUE, SAVE = TRUE, 
                             outputDIR = NULL, fileNAME = NULL, overwrite = TRUE) {
  require(sf); require(dplyr)
  
  ## checks
  if(any(class(sf.obj) != c("sf", "data.frame"))) stop ("sf.obj must be an sf object")
  if(is.null(buff.dist)) stop("Define a buffer distance")
  if(is.null(fireNAMES)) stop("Provide the name of the fire ID variable")
  if(is.null(fireVARS)) stop("Provide a vector of fire variables of interest (indices or variable names)")
  if(SAVE & is.null(outputDIR)) stop("SAVE is TRUE, but output folder is not defined")
  if(SAVE & is.null(fileNAME)) stop("SAVE is TRUE, but file name prefix is not defined")
  
  ## DEFINE PROJECTION AND RE-PROJECT IF NEED BE
  crsProj <- if(is.null(crsProj)) st_crs(sf.obj) else CRS(crsProj)
  
  if(crsProj != st_crs(sf.obj)) {
    warning("Reprojecting sf.obj to selected projection")
    sf.obj <- st_transform(sf.obj, crs = crsProj)
  }
  
  fire.ls <- unique(eval(parse(text = paste0("sf.obj$", fireNAMES))))
  ## CHECK FIRE VARIABLES
  if(is.numeric(fireVARS)) fireVARS <- names(sf.obj)[fireVARS]
  cat(paste0("Using the following fire variables: ", paste0(fireVARS, collapse =", ")), "\n")
  
  fireEvent.ls <- lapply(fire.ls, FUN = function(fire) {
    print(fire)
    browser()
    firePolys <- eval(parse(text = paste0("sf.obj$", fireNAMES))) == fire
    sf.fire <- sf.obj[firePolys, c(fireNAMES, fireVARS)]
    
    ## CALCULATE FIRE AND EVENT PERIMETERS
    firePerim <- st_union(sf.fire$geometry)
    outerBuff <- st_buffer(firePerim, dist = buff.dist)
    eventPerim <- st_buffer(outerBuff, dist = -buff.dist)
    
    ## REMOVE INNER MATRIX HOLES FROM EVENT AND ORIGINAL FIRE PERIMETER
    ## (i.e. unburnt patches surrounded by fire)
    if(class(eventPerim)[1] == "sfc_MULTIPOLYGON") {
      noHolesEventPerim <- st_sfc(st_multipolygon(lapply(eventPerim[[1]], function(x) x[1])), crs = crsProj)
      noHolesEventPerim <- st_union(noHolesEventPerim)  ## union to account for disturbed patches inside inner matrix, nested within larger disturbed patches
    } else {
      noHolesEventPerim <- st_sfc(st_multipolygon(lapply(eventPerim[1], function(x) x[1])), crs = crsProj)
      noHolesEventPerim <- st_union(noHolesEventPerim) ## union to account for disturbed patches inside inner matrix, nested within larger disturbed patches
    }
    
    if(class(firePerim)[1] == "sfc_MULTIPOLYGON") {
      noHolesFirePerim <- st_sfc(st_multipolygon(lapply(firePerim[[1]], function(x) x[1])), crs = crsProj)
      noHolesFirePerim <- st_union(noHolesFirePerim)  
    } else {
      noHolesFirePerim <- st_sfc(st_multipolygon(lapply(firePerim[1], function(x) x[1])), crs = crsProj)
      noHolesFirePerim <- st_union(noHolesFirePerim)
    }
    
    ## EXTRACT INNER MATRIX HOLES PRODUCED BY BUFFERING - these will determine what the inner matrix remnants are
    bufferedHoles <- st_difference(noHolesEventPerim, eventPerim)
    
    ## EXTRACT ALL REMNANTS PRODUCED BY BUFFERING (except holes)
    allRemnants <- st_difference(eventPerim, firePerim)
    
    ## CALCULATE INTERSECTIONS (TOUCH) BETWEEN HOLES AND REMNANTS
    if(length(bufferedHoles) > 0) {
      if(class(allRemnants)[1] == "sfc_MULTIPOLYGON") {
        remnHoleInters  <- sapply(allRemnants[[1]], FUN = function(sfgpoly) {
          sfcpoly <- st_sfc(list = st_polygon(sfgpoly[1]), crs = crsProj)
          st_intersects(sfcpoly, bufferedHoles, sparse = FALSE) 
        })
      } else {
        remnHoleInters  <- sapply(allRemnants[1], FUN = function(sfgpoly) {
          sfcpoly <- st_sfc(list = st_polygon(sfgpoly[1]), crs = crsProj)
          st_intersects(sfcpoly, bufferedHoles, sparse = FALSE) 
        })
      }
    } else {
      remnHoleInters <- if(class(allRemnants)[1] == "sfc_MULTIPOLYGON") {
        rep(FALSE, length(allRemnants[[1]]))
        } else {
        rep(FALSE, length(allRemnants[1]))
      }
    }
    
    ## DEFINE INNER MATRIX REMNANTS
    if(class(allRemnants)[1] == "sfc_MULTIPOLYGON") {
      inMatrixRemn <- lapply(allRemnants[[1]][remnHoleInters], st_polygon) %>% st_sfc(., crs = crsProj)
      inMatrixRemn <- st_union(inMatrixRemn)
    } else {
      inMatrixRemn <- lapply(allRemnants[1][remnHoleInters], st_polygon) %>% st_sfc(., crs = crsProj)  
      inMatrixRemn <- st_union(inMatrixRemn)
    }
    
    ## DEFINE OUTER MATRIX REMNANTS
    outMatrixRemn <- st_difference(allRemnants, inMatrixRemn)
    
    ## MERGE HOLES AND INNER MATRIX
    inMatrixRemn2 <- c(bufferedHoles, inMatrixRemn)
    
    ## CONVERT TO MULTIPOLYGONS AND SIMPLE FEATURE COLLECTION, ADDING FIRE NAME
    firePerim <- st_sfc(firePerim) 
    firePerim <- st_sf(geometry = firePerim)   ## first "combine" list of polygons into a multipolygon, which is then converted to a Simple Features object
    firePerim$PatchType <- "disturbedPatch"
    ## add fire details
    firePerim <- st_intersection(firePerim, sf.fire)
    
    eventPerim <- st_sfc(eventPerim)
    eventPerim <- st_sf(geometry = eventPerim)  
    eventPerim$PatchType <- "eventPerim"
    ## add fire details
    eventPerim <- st_intersection(eventPerim, sf.fire)
    
    outMatrixRemn <- st_sfc(outMatrixRemn)
    outMatrixRemn <- st_sf(geometry = outMatrixRemn) 
    outMatrixRemn$PatchType <- "outMatrixRemn"
    
    inMatrixRemn2 <- st_sfc(inMatrixRemn2)
    inMatrixRemn2 <- st_sf(geometry = inMatrixRemn2) 
    inMatrixRemn2$PatchType <- "inMatrixRemn"
    
    temp.df <- firePerim[, !names(firePerim) %in% names(outMatrixRemn), drop = TRUE]
    temp.df$SEV_CLAS <- NA
    temp.df <- temp.df[!duplicated(temp.df),]
    
    outMatrixRemn <- merge(outMatrixRemn, temp.df)
    inMatrixRemn2 <- merge(inMatrixRemn2, temp.df)
    
    ## COMBINE INTO ONE OBJECT
    fireEvent <- rbind(firePerim, eventPerim, outMatrixRemn, inMatrixRemn2)
    
    return(fireEvent)
  })
  
  fireEvent.all <- do.call(rbind, fireEvent.ls)
  
  if(PLOT) plot(fireEvent.all, key.pos = 1)
  
  ## SAVE AS SHAPEFILE
  if(!dir.exists(outputDIR)) dir.create(outputDIR, recursive = TRUE)
  if(SAVE) st_write(fireEvent.all, file.path(outputDIR, paste0(fileNAME, ".shp")), delete_layer = overwrite)
  
  return(fireEvent.all)
}
