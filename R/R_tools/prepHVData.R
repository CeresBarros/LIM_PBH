## LOAD HV DATA FOR ANALYSES -----------------------

allFiles <- list.files(HVoutputPath, "Intersection.*.rds", full.names = TRUE)
if (mergeDMCPSME) {
  allFiles <- grep("_DMCPSME_|_PSME_|_dryPSME_", allFiles, value = TRUE, invert = TRUE) ## remove HV for vegTypes that were merged
  allFiles <- c(allFiles, list.files(HVoutputPathMergedVegType, "Intersection.*.rds", full.names = TRUE)) ## add HV for merged vegType
}

if (mergePSME) {
  allFiles <- grep("_PSME_|_dryPSME_", allFiles, value = TRUE, invert = TRUE) ## remove HV for vegTypes that were merged
  allFiles <- c(allFiles, list.files(HVoutputPathMergedVegType, "Intersection.*.rds", full.names = TRUE)) ## add HV for merged vegType
}


fireHVData <- loadHVResultsFromRDS("fireHVs", allFiles)
if ("HVid" %in% names(fireHVData)) {
  fireHVData[, scenario := HVid]   ## no HVid with intersecion results
}
fireHVData[, year := as.integer(max(yearSubset))]   ## fire HV are from last year, but integrate the whole simulation period
## drop unique components
set(fireHVData, NULL, grep("Unique", names(fireHVData)), NULL)

## add comparison type
comp <- sub(".*_", "", grep("Volume", names(fireHVData), value = TRUE))
comp <- paste(comp, collapse = "_")
fireHVData[, compare := comp]
setnames(fireHVData, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(fireHVData)))

## check reps/years/scenario per veg type
if (getOption("LandR.assertions")) {
  temp <- split(fireHVData[, .(scenario, year, rep, repHV, vegType)], by = "vegType", keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) paste0(x$scenario, x$year, x$rep))
  temp <- lapply(temp, FUN = unique)
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))
  test <- sapply(test, length)

  if (any(test)) {
    stop("fireHVData has different combinations of scenario, year, rep across vegTypes")
  }
}

## get between and within year (between scenarios) comparisons separately
withinYearFiles <- grep("yr", allFiles, value = TRUE)
vegHVDataWYrComparisons <- loadHVResultsFromRDS("vegHVs", withinYearFiles)
if ("HVid" %in% names(vegHVDataWYrComparisons)) {
  vegHVDataWYrComparisons[, tempHVid := as.numeric(HVid)]
  vegHVDataWYrComparisons[is.na(tempHVid), scenario := HVid]
  vegHVDataWYrComparisons[!is.na(tempHVid), year := HVid]
  vegHVDataWYrComparisons[, tempHVid := NULL]
}
## drop unique components
set(vegHVDataWYrComparisons, NULL, grep("Unique", names(vegHVDataWYrComparisons)), NULL)

## add comparison type
comp <- sub(".*_", "", grep("Volume", names(vegHVDataWYrComparisons), value = TRUE))
comp <- paste(comp, collapse = "_")
vegHVDataWYrComparisons[, compare := comp]

## break into 2 tables, remove empty volume columns,
## make column of comparison ID, change names and re-rbind
tempData  <- vegHVDataWYrComparisons[is.na(Volume_HV1_PM),]
tempData2  <- vegHVDataWYrComparisons[!is.na(Volume_HV1_PM),]

cols <- grep("Volume", names(tempData), value = TRUE)
set(tempData, NULL, cols[which(is.na(colSums(tempData[, ..cols])))], NULL)
set(tempData2, NULL, cols[which(is.na(colSums(tempData2[, ..cols])))], NULL)

setnames(tempData, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData)))
setnames(tempData2, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData2)))

vegHVDataWYrComparisons <- rbind(tempData, tempData2, use.names = TRUE)

## between year comparisons data
betweenYearFiles <- grep("yr", allFiles, value = TRUE, invert = TRUE)
vegHVDataBYrComparisons <- loadHVResultsFromRDS("vegHVs", betweenYearFiles)
if ("HVid" %in% names(vegHVDataBYrComparisons)) {
  vegHVDataBYrComparisons[, tempHVid := as.numeric(HVid)]
  vegHVDataBYrComparisons[is.na(tempHVid), scenario := HVid]
  vegHVDataBYrComparisons[!is.na(tempHVid), year := HVid]
  vegHVDataBYrComparisons[, tempHVid := NULL]
}
## drop unique components
set(vegHVDataBYrComparisons, NULL, grep("Unique", names(vegHVDataBYrComparisons)), NULL)

## sometimes HV1 is 2011, others is 2111 (probably because data wasn't sorted) -
## make sure there aren't duplicate comparisons
if (getOption("LandR.assertions")) {
  test <- unique(vegHVDataBYrComparisons[is.na(Volume_HV1_2011), .(scenario, rep, vegType)])
  test2 <- unique(vegHVDataBYrComparisons[!is.na(Volume_HV1_2011), .(scenario, rep, vegType)])
  test[, test := paste(scenario, rep, vegType)]
  test2[, test := paste(scenario, rep, vegType)]
  if (length(intersect(test$test, test2$test)) |
      length(intersect(test2$test, test$test))) {
    stop("There are duplicated HV comparisions between years\n",
         "(per scenario/rep/vegType combination)")
  }

  test <- any(is.na(vegHVDataBYrComparisons[is.na(Volume_HV1_2011), Volume_HV1_2111]))
  test2 <- any(is.na(vegHVDataBYrComparisons[!is.na(Volume_HV1_2011), Volume_HV1_2111]))
  test3 <- any(is.na(vegHVDataBYrComparisons[is.na(Volume_HV2_2011), Volume_HV2_2111]))
  test4 <- any(is.na(vegHVDataBYrComparisons[!is.na(Volume_HV2_2011), Volume_HV2_2111]))

  if ((isTRUE(test) | isFALSE(test2)) |
      (isTRUE(test3) | isFALSE(test4))) {
    stop("There are either duplicated or missing HV comparisions between years\n",
         "(per scenario/rep/vegType combination)")
  }
}

## add comparison type - because sometimes HV1 is 2011, others is 2111
## we need to use `sort(unique())` bellow
comp <- sort(unique(sub(".*_", "", grep("Volume", names(vegHVDataBYrComparisons), value = TRUE))))
comp <- paste(comp, collapse = "_")
vegHVDataBYrComparisons[, compare := comp]

## break into 2 tables, remove empty volume columns,
## make column of comparison ID, change names and re-rbind
tempData  <- vegHVDataBYrComparisons[is.na(Volume_HV1_2011),]
tempData2  <- vegHVDataBYrComparisons[!is.na(Volume_HV1_2011),]

cols <- grep("Volume", names(tempData), value = TRUE)
set(tempData, NULL, cols[which(is.na(colSums(tempData[, ..cols])))], NULL)
set(tempData2, NULL, cols[which(is.na(colSums(tempData2[, ..cols])))], NULL)

setnames(tempData, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData)))
setnames(tempData2, new = sub("Volume_(HV)[0-9]_(.*)", "\\1_\\2", names(tempData2)))

vegHVDataBYrComparisons <- rbind(tempData, tempData2, use.names = TRUE)

vegHVData <- rbind(vegHVDataWYrComparisons, vegHVDataBYrComparisons, use.names = TRUE,
                   fill = TRUE)

## check reps/years/scenario per veg type
if (getOption("LandR.assertions")) {
  temp <- split(vegHVData[, .(scenario, year, rep, vegType)], by = "vegType", keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) paste0(x$scenario, x$year, x$rep))
  temp <- lapply(temp, FUN = unique)
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))
  test <- sapply(test, length)

  if (any(test)) {
    stop("vegHVData has different combinations of scenario, year, rep across vegTypes")
  }
}

## bind the two tables
allHVData <- rbindlist(list("fireHV" = fireHVData, "vegHV" = vegHVData),
                       use.names = TRUE, fill = TRUE, idcol = "HVtype")
allHVData$vegType <- as.factor(allHVData$vegType)

## calculate overlap following Barros et al 2016
allHVData[, overlap := Intersection/Union]

## check for missing data
if (nrow(allHVData[(is.na(HV_noPM)|is.na(HV_PM)) & (is.na(HV_2011)|is.na(HV_2111))])) {
  stop("There seems to be missing data")
}

rm(temp, tempData, tempData2, fireHVData, vegHVData)
gc(reset = TRUE)
