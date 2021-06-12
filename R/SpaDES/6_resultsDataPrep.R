## -----------------------------------
##  DATA PREP FOR ANALYSES OF RESULTS
## -----------------------------------

## GET CAMERON'S AGE DATA AND STAND VEG TYPES -----------------------
ageDataCN <- fread("data/CameronsAgeData/treelist_outputs_for Ceres.csv")
patchVegTypeCN <- fread("data/CameronsAgeData/patch outputs_for Ceres.csv")

ageDataCN <- patchVegTypeCN[ageDataCN, on = .(Patch.ID)]
ageDataCN$Cover.dendro <- sub("Mixedwood", "mixedwood", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("Broadleaf", "broadleaf", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("-", "", ageDataCN$Cover.dendro)
rm(patchVegTypeCN)

## remove a record that seems funky (maybe it's a new cohort?)
ageDataCN <- ageDataCN[Reconstructed.age != 2018]

## calculate no. cohorts:
ageDataCN[, noCohorts := length(unique(Est.bin)) , by = .(Patch.ID)]


## GET OUTPUTS FOLDERS FOR EACH SCENARIO -----------------------------
outputs_PM <- list.dirs(simPaths$outputPath, full.names = TRUE, recursive = TRUE) %>%
  grep("/PM_rep[[:digit:]]*$", ., value = TRUE)
outputs_noPM <- list.dirs(simPaths$outputPath, full.names = TRUE, recursive = TRUE) %>%
  grep("/noPM_rep[[:digit:]]*$", ., value = TRUE)

## GET FILE NAMES
cohortDataFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "cohortData", full.names = TRUE)),
                     sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "cohortData", full.names = TRUE))) %>%
  unique(.)

rstCurrentFiresFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "rstCurrentFires", full.names = TRUE)),
                          sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "rstCurrentFires", full.names = TRUE))) %>%
  unique(.)

pixelGroupMapFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "pixelGroupMap", full.names = TRUE)),
                        sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "pixelGroupMap", full.names = TRUE))) %>%
  unique(.)

vegTypeMapFiles <- c(sapply(outputs_PM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE)),
                     sapply(outputs_noPM, FUN = function(x) list.files(x, pattern = "vegTypeMap", full.names = TRUE))) %>%
  unique(.)

## load rasters as stacks
rstCurrentFiresStk_noPM <- lapply(grep("noPM", rstCurrentFiresFiles, value = TRUE), readRDS) %>%
  stack(.)
rstCurrentFiresStk_PM <- lapply(grep("noPM", rstCurrentFiresFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)
pixelGroupMapStk_noPM <- lapply(grep("noPM", pixelGroupMapFiles, value = TRUE), readRDS) %>%
  stack(.)
pixelGroupMapStk_PM <- lapply(grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)
vegTypeMapStk_noPM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE), readRDS) %>%
  stack(.)
vegTypeMapStk_PM <- lapply(grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE), readRDS) %>%
  stack(.)

names(rstCurrentFiresStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentFiresFiles, value = TRUE))),
                                        sub(".*(rep)([0-9]+)/.*", "\\1\\2",  grep("noPM", rstCurrentFiresFiles, value = TRUE)), sep = "_")

names(rstCurrentFiresStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", rstCurrentFiresFiles, value = TRUE, invert = TRUE))),
                                      sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", rstCurrentFiresFiles, value = TRUE, invert = TRUE)), sep = "_")
names(pixelGroupMapStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE))),
                                      sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", pixelGroupMapFiles, value = TRUE)), sep = "_")

names(pixelGroupMapStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE))),
                                    sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", pixelGroupMapFiles, value = TRUE, invert = TRUE)), sep = "_")

names(vegTypeMapStk_noPM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE))),
                                   sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", vegTypeMapFiles, value = TRUE)), sep = "_")

names(vegTypeMapStk_PM) <- paste(sub(".*year", "year", sub("\\.rds", "", grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE))),
                                 sub(".*(rep)([0-9]+)/.*", "\\1\\2", grep("noPM", vegTypeMapFiles, value = TRUE, invert = TRUE)), sep = "_")



yearSubset <- c(seq(1, 100, 5), 100)

## BUILD TABLES OF RESULTS ------------------------
## pixelCohortData tables
files <- grep("noPM", cohortDataFiles, value = TRUE)
pixelCohortData_noPM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  r <- sub(".*(rep)([0-9]+)/.*", "\\2", ff)
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr, "_rep", r)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  cohortData[, rep := as.integer(r)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_noPM) %>%
  rbindlist(fill = TRUE, l = .)
pixelCohortData_noPM <- pixelCohortData_noPM[year %in% yearSubset]

files <- grep("noPM", cohortDataFiles, value = TRUE, invert = TRUE)
pixelCohortData_PM <- lapply(files, FUN = function(ff, pixelGroupMapStk) {
  cohortData <- readRDS(ff)
  yr <- sub(".*year", "",  sub(".rds", "", ff))
  r <- sub(".*(rep)([0-9]+)/.*", "\\2", ff)
  pixelGroupMap <- pixelGroupMapStk[[paste0("year", yr, "_rep", r)]]
  cohortData <- addPixels2CohortData(cohortData, pixelGroupMap)
  cohortData[, year := as.integer(yr)]
  cohortData[, rep := as.integer(r)]
  return(cohortData)
}, pixelGroupMapStk = pixelGroupMapStk_PM) %>%
  rbindlist(fill = TRUE, l = .)
pixelCohortData_PM <- pixelCohortData_PM[year %in% yearSubset]

## vegTypeData tables
vegTypeSubset <- intersect(names(vegTypeMapStk_noPM), names(pixelGroupMapStk_noPM))
vegTypeData_noPM <- lapply(vegTypeSubset, FUN = function(x) {
  yr <- grep("year", unlist(strsplit(x, split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(x, split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_noPM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
             vegType = vegTypeMapStk_noPM[[x]][],
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
}) %>%
  rbindlist(.)
vegTypeData_noPM <- vegTypeData_noPM[!is.na(pixelGroup)]
vegTypeData_noPM <- vegTypeData_noPM[year %in% yearSubset]

vegTypeSubset <- intersect(names(vegTypeMapStk_PM), names(pixelGroupMapStk_PM))
vegTypeData_PM <- lapply(vegTypeSubset, FUN = function(x) {
  yr <- grep("year", unlist(strsplit(x, split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(x, split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(vegTypeMapStk_PM[[x]])),
             pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
             vegType = vegTypeMapStk_PM[[x]][],
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)))
})  %>%
  rbindlist(.)
vegTypeData_PM <- vegTypeData_PM[!is.na(pixelGroup)]
vegTypeData_PM <- vegTypeData_PM[year %in% yearSubset]

## pixelBurnData tables - all rasters
pixelBurnData_noPM <- lapply(unstack(rstCurrentFiresStk_noPM), FUN = function(ras) {
  yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(ras)),
             burnt = as.integer(!is.na(getValues(ras))),
             fireID = as.integer(getValues(ras)),
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)),
             scenario = "noPM")
}) %>%
  rbindlist(.)

pixelBurnData_PM <- lapply(unstack(rstCurrentFiresStk_PM), FUN = function(ras) {
  yr <- grep("year", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  r <- grep("rep", unlist(strsplit(names(ras), split = "_")), value = TRUE)
  data.table(pixelIndex = seq_len(ncell(ras)),
             burnt = as.integer(!is.na(getValues(ras))),
             fireID = as.integer(getValues(ras)),
             year = as.integer(sub("year", "", yr)),
             rep = as.integer(sub("rep", "", r)),
             scenario = "PM")
})  %>%
  rbindlist(.)

allPixelBurnData <- rbind(pixelBurnData_noPM, pixelBurnData_PM, use.names = TRUE)

## clean-up to save memory - keep vegTypeMapStk_noPM for labels
rm(files, vegTypeSubset, rstCurrentFiresStk_PM, rstCurrentFiresStk_noPM,
   vegTypeMapStk_PM, pixelGroupMapStk_PM, pixelGroupMapStk_noPM,
   pixelBurnData_PM, pixelBurnData_noPM)
amc::.gc()

## join tables, add scenario col and rbind.
## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
## pixelBurntData can also have a different number of years if the saving frequency differs
## so join by keeping all pixels, calculate fire properties per pixel, then subset to
## cohort data years.
pixelCohortData_noPM <- merge(vegTypeData_noPM, pixelCohortData_noPM,
                              by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE)
pixelCohortData_noPM[, scenario := "noPM"]
amc::.gc()

pixelCohortData_PM <- merge(vegTypeData_PM, pixelCohortData_PM,
                            by = c("pixelIndex", "pixelGroup", "year", "rep"), all = TRUE)
pixelCohortData_PM[, scenario := "PM"]
amc::.gc()

allPixelCohortData <- rbind(pixelCohortData_noPM, pixelCohortData_PM, use.names = TRUE)
if (exists("allPixelCohortData")) rm(list = grep("^(pixelCohort|vegType)Data", ls(), value = TRUE))
amc::.gc()

## ECOLOGICAL ZONATION -----------------------------
preSimList <- loadSimList(file.path(simPaths$outputPath, "noPM", "LIM_simInit_noPM"))
ecoregionLayerRas <- rasterize(preSimList$ecoregionLayer, preSimList$rasterToMatch, field = "ecozoneCode")
ecoregionLayerDT <- data.table(ecozoneCode = getValues(ecoregionLayerRas),
                               pixelIndex = seq_len(ncell(ecoregionLayerRas)))
ecoregionLayerDT <- ecoregionLayerDT[!is.na(ecozoneCode)]

ecoregionLayerLabels <- data.table(ecozoneCode = preSimList$ecoregionLayer$ecozoneCode,
                                   ecozoneName = paste(preSimList$ecoregionLayer$NRNAME,
                                                       preSimList$ecoregionLayer$NSRNAME, sep = " - ")) %>%
  unique(.)

ecoregionLayerDT <- ecoregionLayerLabels[ecoregionLayerDT, on = .(ecozoneCode)]

allPixelCohortData <- ecoregionLayerDT[allPixelCohortData, on = .(pixelIndex)]
amc::.gc()


## NO. FIRES ---------------------------------------
## how many times did each pixel burn?
## calculate total no. fires per pixel/scenario/rep
## calculate fire size in pixels per fireID/scenario/rep
## calculate fire frequency as the mean fire-intervals per pixel (see Steel et al 2021 for limitations and details)

allPixelBurnData <- allPixelBurnData[fireID != "NA"]
allPixelBurnData[, noFires := sum(burnt), by = .(scenario, rep, pixelIndex)]
allPixelBurnData[, fireSize := length(pixelIndex), by = .(scenario, rep, year, fireID)]

## fire frequency
fireFreqDT <- allPixelBurnData[, list(year = c(year, 100),     ## year one is dropped here
                                      fireInt = diff(c(1, year, 100))),   ## interval calculated between fire years, and start and end years
                               by = .(scenario, rep, pixelIndex)]
## because we forced a start and end year, intervals of 0 for the 100th year mean that there was a fire at year 100
## this doesn't apply in the same way to fires at year 1, these should have a return interval of 0 (because we added the first year)
## if only one fire occurred and it was at year 100, then the correct interval is 99 (100-1)
fireFreqDT <- fireFreqDT[!(fireInt == 0 & year == 100)]
fireFreqDT[, fireFreq := mean(fireInt), by = .(scenario, rep, pixelIndex)]

if (any(is.na(fireFreqDT$fireFreq))) stop("NA fire intervals")

## join DTs
allPixelBurnData <- fireFreqDT[allPixelBurnData, on = .(scenario, rep, year, pixelIndex)]
allPixelBurnData[, burnt := NULL] ## no longer necessary

test <- sapply(split(allPixelBurnData, by = c("scenario", "rep", "year")), FUN = function(x){
  any(duplicated(x[, pixelIndex]))
})
if (any(test))
  stop("Each pixel should only have one record of no. fires per scenario")

rm(fireFreqDT, test)
amc::.gc()

## add noFires to cohortData - no year info, because its the total across the simulation
cols <- c("scenario", "rep", "pixelIndex", "noFires")
allPixelCohortData <- unique(allPixelBurnData[, ..cols])[allPixelCohortData,
                                                         on = .(scenario, rep, pixelIndex)]

cols <- c("scenario", "rep", "pixelIndex")

test <- allPixelBurnData[allPixelCohortData[is.na(noFires), ..cols], on = cols, nomatch = 0]
test <- dim(test)[1]
if (test) stop("There shouldn't be any NAs in noFires except in scenario/rep/pixelIndex\n",
               "combos that did not have any fires during the simulation")

allPixelCohortData[is.na(noFires), noFires := 0]
allPixelCohortData <- allPixelCohortData[!is.na(pixelGroup),]
amc::.gc()

## ADD MISSING SPECIES IN YEAR/SCENARIO/PIXEL COMBINATION
## cohortData doens't track absent cohorts, so they need to ba added back
## for now pixels from the 0s pixelGroup  have one entry with NAs for speciesCode
## they will be ignored for now and removed later, after adding one species entry for each of these pixels.
## for reporting consistency add to show losses in B
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
test <- length(unique(allPixelCohortData[, length(unique(pixelIndex)), by = .(scenario, rep, year)]$V1)) == 1
test2 <- length(unique(allPixelCohortData[, length(unique(speciesCode)), by = .(scenario, rep, year, pixelIndex)]$V1)) == 1
test3 <- any(is.na(allPixelCohortData$speciesCode))

if (isFALSE(test))
  stop("No. pixels should be the same across years, for a given scenario/rep")
if (isFALSE(test2))
  stop("No. species per pixel should be the same across pixels, for a given scenario/rep/year")
if (test3)
  stop("There are NA speciesCodes")



## add ecoregion group/ecozone code/name where they're missing
## add vegType where it's missing, but it's a pixel with some veg
allPixelCohortData[, `:=`(vegType = max(vegType, na.rm = TRUE)),
                   by = .(scenario, rep, year, pixelGroup)]

allPixelCohortData[, `:=`(ecoregionGroup = unique(na.omit(ecoregionGroup)),
                          ecozoneCode = unique(na.omit(ecozoneCode)),
                          ecozoneName = unique(na.omit(ecozoneName))),
                   by = .(pixelIndex)]
amc::.gc()

## add noFires where it's missing
allPixelCohortData[, noFires := max(noFires, na.rm = TRUE),
                   by = .(scenario, rep, pixelIndex)]
amc::.gc()

## replace NAs of cohortData by 0s
cols <- c("age", "B", "mortality", "aNPPAct", "vegType", "noFires")
allPixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
amc::.gc()

## add presence/absence of fire across simulation per pixel/scenario
allPixelCohortData[, firePresAbs := as.integer(any(noFires > 0)), by = .(scenario, rep, pixelIndex)]
amc::.gc()


## USING CAMERON'S CLASSIFICATION/SUMMARY ---------------------
## Cameron uses relative basal area to classify stand structure, we can use relative Biomass.
## we subset to the montane ecological zone, from where Cameron's data comes from
allPixelCohortDataMnt <- allPixelCohortData[grep("Montane", ecozoneName)]
allPixelCohortDataMnt[, sumB := sum(B), by = .(scenario, year, rep, pixelGroup)]
allPixelCohortDataMnt[, relB := sum(B)/sumB, by = .(scenario, rep, year, pixelGroup, speciesCode)]
allPixelCohortDataMnt[is.na(relB) & sumB == 0, relB := 0]
amc::.gc()

if (any(is.na(allPixelCohortDataMnt$relB)))
  stop("Missing values in relative biomass")

## subset to a smaller DT and join Cameron's species names
vegTypesCN <- unique(allPixelCohortDataMnt[B > 0, .(scenario, rep, year, pixelGroup, speciesCode, relB)])
vegTypesCN <- unique(na.omit(preSimList$sppEquiv[, .(Cameron, LIM)]))[vegTypesCN, on = "LIM==speciesCode",
                                                                      allow.cartesian = TRUE]
setnames(vegTypesCN, "LIM", "speciesCode")

parallelFUN <- function(DT) {
  set.seed(123)
  tempArg <- sample(1:nrow(DT), 100, replace = FALSE)
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
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = 8)
} else {
  plan("multicore", workers = 8)
}

vegTypesCN <- future_lapply(split(vegTypesCN, by = c("scenario", "rep")),
                            FUN = parallelFUN)
future:::ClusterRegistry("stop")
amc::.gc()

vegTypesCN <- rbindlist(vegTypesCN, use.names = TRUE)

## test:
# vegTypesCN <- lapply(unique(showCache(cPath,
#                                       userTags = c("convertToCNVegType", "reportResults"),
#                                       after = "2021-05-05")$cacheId),
#                      FUN = function(x) {loadFromCache(cPath, cacheId = x)}) %>%
#   rbindlist(.)

## add Cameron's veg types and get rid of useless columns
cols <- c("scenario", "rep", "year", "pixelGroup", "vegTypeCN")
cols2 <- c("scenario", "rep", "year", "pixelGroup")
allPixelCohortDataMnt <- tryCatch(unique(vegTypesCN[, ..cols])[allPixelCohortDataMnt,
                                                               on = cols2],
                                  error = allPixelCohortDataMnt)
if (!"vegTypeCN" %in% names(allPixelCohortDataMnt))
  stop("Joining Cameron's veg types didn't work") else
    rm(vegTypesCN)

if (any(is.na(allPixelCohortDataMnt$vegTypeCN) & allPixelCohortDataMnt$B > 0))
  stop("Some pixels with biomass were not assigned a vegTypeCN")

allPixelCohortDataMnt[is.na(vegTypeCN), vegTypeCN := "No veg."]
allPixelCohortDataMnt[, `:=`(sumB = NULL,
                             relB = NULL,
                             vegType = NULL)]
## make "No veg." the last factor
levs <- c(sort(grep("No veg.", unique(allPixelCohortDataMnt$vegTypeCN), value = TRUE, invert = TRUE)),
          "No veg.")
allPixelCohortDataMnt[, vegTypeCN := factor(vegTypeCN, levels = levs)]
amc::.gc()
