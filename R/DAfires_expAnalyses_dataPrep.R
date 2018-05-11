## ---------------------------------
## DAVID ANDISON FIRE DATA
##
## Dataprep for Exploratory analysis
## ---------------------------------

## this script should be sourced

rm(list = ls()); gc(reset = TRUE)

## using - as of April 19th, 2018
# loading reproducible     0.1.4.9015
# loading quickPlot        0.1.3.9002
# loading SpaDES.core      0.1.1.9009
# loading SpaDES.tools     0.1.1.9005
# loading SpaDES.addins    0.1.1
# devtools::install_github("PredictiveEcology/reproducible@development")
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")

## requires
library(SpaDES); library(sf); library(ggplot2); library(data.table); library(dplyr)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## LOAD DATA ---------------------------------------

## GET FIRE DATA
## THIS IS NOT WORKING:
# dataDIR <- "data/fires_Dave" 
# DA_fires <- prepInputs(targetFile = file.path(dataDIR, "albertafires1_postfire.shp"),
#                        url = "https://drive.google.com/file/d/1wGCqI_X1t-PDM5eO6JWQW9hpNv_8zlum/view?usp=sharing",
#                        archive = "Projected_renamed.zip",
#                        alsoExtract = file.path(dataDIR, c("albertafires1_postfire.", "albertafires2_postfire.", 
#                                                           "albertafires1_pretfire.", "albertafires2_prefire.",
#                                                           "saskatchewanfires_postfire.","saskatchewanfires_prefire.")),
#                        fun = "st_read", pkg = "sf",
#                        destinationPath = dataDIR)

files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire")
folder = "~/../OneDrive/Documents/LandscapesInMotion/data/fires_Dave/Projected_renamed"

for(x in files) {
  eval(parse(text = paste0(
    x, " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
    )))
}

head(albertafires1_postfire)
head(albertafires2_postfire)
head(saskatchewanfires_postfire)

## Saskatchewan attributes table has some repeated columns, check if they are duplicated before deleting
if(all(saskatchewanfires_postfire$FIRE_NAME_ == saskatchewanfires_postfire$FIRE_NAME)) {
  saskatchewanfires_postfire$FIRE_NAME_ <- NULL
}
if(all(saskatchewanfires_postfire$AREA_ == saskatchewanfires_postfire$AREA)) {
  saskatchewanfires_postfire$AREA_ <- NULL
}
if(all(saskatchewanfires_postfire$PERIMETER_ == saskatchewanfires_postfire$PERIMETER)) {
  saskatchewanfires_postfire$PERIMETER_ <- NULL
}

## Uniformize column names
names(albertafires2_postfire)[grep("YEAR", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("EAR", names(saskatchewanfires_postfire))] = 
  grep("YEAR", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("NAME", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("NAME", names(saskatchewanfires_postfire))] = 
  grep("NAME", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("FIRE_CODE", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("FIRE_NUM", names(albertafires1_postfire))] = 
  grep("FIRE_ID", names(saskatchewanfires_postfire), value = TRUE)

names(albertafires2_postfire)[grep("ELEMENT", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("SURVIVAL", names(albertafires1_postfire))] = 
  names(saskatchewanfires_postfire)[grep("CLASS", names(saskatchewanfires_postfire))] = "SEV_CLAS"

## PRE FIRE DATA
files = c("albertafires1_prefire", "albertafires2_prefire", "saskatchewanfires_prefire")
folder = "~/../OneDrive/Documents/LandscapesInMotion/data/fires_Dave/Projected_renamed"

for(x in files) {
  eval(parse(text = paste0(
    x, " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
  )))
}

## WATER DATA
# files = c("water-abta", "water-sask")
# folder = "data/fires_Dave/Projected_renamed"
# 
# for(x in files) {
#   eval(parse(text = paste0(
#     sub("-", "_", x), " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
#   )))
# }

## DEFINE FIRE EVENTS ----------------------------------------

# firesABSK <- Cache(loadBindSpatialObjs, 
#                    files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire"),
#                    folder = "data/fires_Dave/Projected_renamed", 
#                    cacheRepo = "analyses/cache", userTags = "fireData")

## Use Alberta 1 post fire data only for now, as severity classes on other datasets and not yet comparable.
rm(albertafires2_postfire, albertafires2_prefire, saskatchewanfires_postfire, saskatchewanfires_prefire)
gc(reset = TRUE)
AB1_fireEvents <- Cache(defineFireEvents, 
                    sf.obj = albertafires1_postfire, fireNAMES = "FIRE_NAME",
                    # fireVARS = c("FIRE_ID", "FIRE_YEAR", "SEV_CLAS"),   ## this makes the output object huge
                    buff.dist = 200L, 
                    PLOT = FALSE, SAVE = FALSE, outputDIR = "analyses/FireEvents", 
                    fileNAME = "Andison_AB1_fireEvents", overwrite = TRUE,
                    cacheRepo = "analyses/cache", userTags = "dataTreat_fireEvents",
                    omitArgs = c("PLOT", "SAVE", "outputDIR", "fileNAME", "overwrite"))

## ADD SEVERITY
## (doing it in defineFireEvents seems to produce an overly large polygon)
AB1_distPatchSev <- Cache(st_intersection, 
                          x = AB1_fireEvents[AB1_fireEvents$PatchType == "disturbedPatch",], 
                          y = albertafires1_postfire[, "SEV_CLAS"], userTags = "dataTreat_fireEvents",
                          cacheRepo = "analyses/cache")

## make NA columns for binding
AB1_fireEvents[, setdiff(names(AB1_distPatchSev), names(AB1_fireEvents))] <- NA
AB1_fireEventsSev <- rbind(AB1_fireEvents[AB1_fireEvents$PatchType != "disturbedPatch",], AB1_distPatchSev)

## JOIN WATER, VEGETATION AND FIRE EVENTS --------------------
AB1_vegFireEvents <- Cache(st_intersection, 
                           x = albertafires1_prefire, y = AB1_fireEventsSev,
                           userTags = "dataTreat_fireEvents_wVeg",
                           cacheRepo = "analyses/cache")

## save - not working, R thinks this is sfc instead of sf.
# st_write(st_as_sf(AB1_vegFireEvents), 
#          dsn = "analyses/FireEvents/Andison_AB1_fireEventsVegetation.shp",
#          delete_layer = TRUE)  ##  not working
# raster::shapefile(as_Spatial(st_as_sf(AB1_vegFireEvents)),
#                   filename = "analyses/FireEvents/Andison_AB1_fireEventsVegetation.shp",
#                   overwrite = TRUE)

## not sure what to do about water yet... 
## perhaps just remove these areas with st difference before intersecting with veg?
# AB1_watVegFireInters <- Cache(st_intersection,
#                               x = AB1_vegFireEvents, y = water_abta[, "FEATURE_TY"],
#                               userTags = "dataTreat_fireEvents_wVeg_water", cacheRepo = "analyses/cache")
# names(AB1_watVegFireEvents)[which(names(AB1_watVegFireEvents) == "FEATURE_TY")] <- "WATER_TY"

## MAKE DATATABLE ----
## extract DATATABLE only
AB1_vegFireEvents.dt <- as.data.table(st_set_geometry(AB1_vegFireEvents, NULL))

## remove columns with NAs only
NAcols <- sapply(AB1_vegFireEvents.dt, FUN = function(x) {
  return(any(sum(is.na(x)) == length(x)))
})

## remove duplicates
AB1_vegFireEvents.dt <- AB1_vegFireEvents.dt[!duplicated(AB1_vegFireEvents.dt),]

## MAKE SEVERITY CONTINUOUS --------------

## alberta 1 classes are in survival %, invert scale
AB1_vegFireEvents.dt[SEV_CLAS == "100%", SEV_CONT:= 0]
AB1_vegFireEvents.dt[SEV_CLAS == "75-99%", SEV_CONT:= median(c(1,25))]
AB1_vegFireEvents.dt[SEV_CLAS == "50-74%", SEV_CONT:= median(c(26,50))]
AB1_vegFireEvents.dt[SEV_CLAS == "25-49%", SEV_CONT:= median(c(51,75))]
AB1_vegFireEvents.dt[SEV_CLAS == "6-24%", SEV_CONT:= median(c(76,94))]
AB1_vegFireEvents.dt[SEV_CLAS == "0-5%", SEV_CONT:= median(c(95,100))]

## attribute 0 severity to remnants
AB1_vegFireEvents.dt[PatchType %in% c("outMatrixRemn", "inMatrixRemn"), SEV_CONT:= 0]

## MELT DATA -------
id.vars <- grep("SP", names(AB1_vegFireEvents.dt), value = TRUE, invert = TRUE)

dt.ls <- lapply(grep("SP.$", names(AB1_vegFireEvents.dt), value = TRUE), FUN = function(x) {
  measure.vars <- grep(paste0("^", x), names(AB1_vegFireEvents.dt), value = TRUE)
  temp <- cbind(AB1_vegFireEvents.dt[, id.vars, with = FALSE], AB1_vegFireEvents.dt[, measure.vars, with = FALSE])
  temp <- temp[!duplicated(temp),]
  temp <- temp[!is.na(get(x))]
  
  if(dim(temp)[1] != 0) {
    form <- paste("... ~", measure.vars[1])
    temp2 <- dcast(temp, formula = as.formula(form), value.var = measure.vars[2], fill = 0)
    old.names =  setdiff(names(temp2), names(temp))
    new.names = paste0(measure.vars[1], "_", setdiff(names(temp2), names(temp)))
    names(temp2)[names(temp2) %in% old.names] <- new.names
    return(temp2)
  }
})

## remove null elements
dt.ls <- dt.ls[!sapply(dt.ls, is.null)]
## merge data.tables with Reduce 
AB1_vegFireEvents.mdt <- as.data.table(Reduce(f = function(x,y) merge(x, y, by = id.vars, all = TRUE), dt.ls))

## Replace NAs in spp percentages by 0s (true NAs are not present anymore form subsetting above)
SP.vars <- grep("SP", names(AB1_vegFireEvents.mdt), value = TRUE)
for (j in SP.vars) {
  AB1_vegFireEvents.mdt[is.na(get(j)), (j):=0]
}


## SUM SPP COVER ACROSS DOMINANCE STATUS -------
## cover classes were divieded by 10 so that 100% is = 1 (not 10)

## Canopy
SP.vars <- grep("^SP", names(AB1_vegFireEvents.mdt), value = TRUE)
SP.names <- unique(sub(".*_", "", SP.vars))

SP.cover <- do.call(cbind, lapply(SP.names, FUN = function(x) {
  cols <-  grep(x, SP.vars, value = TRUE)
  temp <- data.table(x = rowSums(AB1_vegFireEvents.mdt[, cols, with = FALSE])/10)
  names(temp) = paste0("SP_", x)
  return(temp)
}))

#Understory
USP.vars <- grep("USP", names(AB1_vegFireEvents.mdt), value = TRUE)
USP.names <- unique(sub(".*_", "", USP.vars))

USP.cover <- do.call(cbind, lapply(USP.names, FUN = function(x) {
  cols <-  grep(x, USP.vars, value = TRUE)
  temp <- data.table(x = rowSums(AB1_vegFireEvents.mdt[, cols, with = FALSE])/10)
  names(temp) = paste0("USP_", x)
  return(temp)
}))

## Merge and add remaining variables
AB1_vegFireEvents.mdt <- cbind(AB1_vegFireEvents.mdt[, id.vars, with = FALSE], cbind(SP.cover, USP.cover))

## CLEAN WS 
rm(dt.ls, id.vars, j, NAcols, SP.vars, SP.names,
   USP.vars, USP.names, SP.cover, USP.cover)

## CALCULATE RELATIVE OCCURRENCES OF VEGETATION ATTRIBUTES PER PATCH/FIRE EVENT -----------
# lapply(grep("SP", names(AB1_vegFireEvents.dt), value = TRUE), FUN = function(x) {
#   browser()
#   temp.dt <- AB1_vegFireEvents.dt %>%
#     group_by_at(c("FIRE_NAME", "PatchType", x)) %>% 
#     summarise(Count = n()) %>%
#     group_by(., FIRE_NAME, PatchType) %>%
#     mutate(.data = ., totalCount = sum(Count)) %>%
#     mutate(.data = ., Freq = Count/totalCount)
#   
#   temp.dt <- melt(temp.dt, id.vars = setdiff(names(temp.dt), c(x)), 
#                   variable.name = "Attribute")
#   temp.dt
#   
# })

