## VEGETATION DATA FOR HVs -----------------------

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't actually change anything because we already subset the pixels with veg
## in the montane belt (regardless of fire)

vegDataForHVs <- allPixelCohortDataMnt[year %in% c(start(preSimList), end(preSimList))]

if (getOption("LandR.assertions")) {
  pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
  temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
  setkey(temp, scenario, rep, pixelIndex)
  setkey(vegDataForHVs, scenario, rep, pixelIndex)

  if (isFALSE(identical(temp, vegDataForHVs))) {
    stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
         "Combinations as allPixelCohortDataMnt")
  }

  temp <- vegDataForHVs[, length(unique(pixelIndex)), by = .(scenario, rep, year)]
  if (unique(temp$V1) > 1) {
    stop("There should be the same number of pixels every year.")
  }

  temp <- split(vegDataForHVs[, .(pixelIndex, scenario, rep, year)],
                by = c("scenario", "rep", "year"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) unique(x[["pixelIndex"]]))
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))

  test <- sapply(test, length)
  if (any(test))
    stop("Different pixelIndex between scenario/rep/year combinations")

  temp <- split(vegDataForHVs[year == start(preSimList), .(vegTypeCN, pixelIndex, scenario, rep)],
                by = c("scenario", "rep"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) setkey(x, vegTypeCN, pixelIndex))
  temppix <- lapply(temp, FUN = function(x) x[["pixelIndex"]])
  tempveg <- lapply(temp, FUN = function(x) x[["vegTypeCN"]])

  test <- lapply(1:length(temppix), function(n) setdiff(temppix[[n]], unlist(temppix[-n])))
  test2 <- lapply(1:length(tempveg), function(n) setdiff(tempveg[[n]], unlist(tempveg[-n])))
  test <- sapply(test, length)
  test2 <- sapply(test2, length)

  if (any(test) | any(test2))
    stop("Difference pixelIndex/vegTypeCN combinations between scenario/reps in the first year")


  ## checks at landscape scale:
  pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
  pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
  temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
  setkey(temp, scenario, rep, pixelIndex)
  setkey(vegDataForHVs, scenario, rep, pixelIndex)

  if (isFALSE(identical(temp, vegDataForHVs))) {
    stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
         "Combinations as allPixelCohortDataMnt")
  }

  temp <- vegDataForHVs[, length(unique(pixelIndex)), by = .(scenario, rep, year)]
  if (unique(temp$V1) > 1) {
    stop("There should be the same number of pixels every year.")
  }

  temp <- split(vegDataForHVs[, .(pixelIndex, scenario, rep, year)],
                by = c("scenario", "rep", "year"), keep.by = FALSE)
  temp <- lapply(temp, FUN = function(x) unique(x[["pixelIndex"]]))
  test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))

  test <- sapply(test, length)
  if (any(test))
    stop("Different pixelIndex between scenario/rep/year combinations")
}

## prep data for hypervolumes
## calculate stand-age as the mean biomass-weighted age
vegDataForHVs[, `:=`(meanStandAge = mean(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                      sum((B/100), na.rm = TRUE))),
                     sdStandAge = sd(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                  sum((B/100), na.rm = TRUE)))),
              by = c("scenario", "rep", "year", "pixelIndex", "vegTypeCN")]
vegDataForHVs[is.na(meanStandAge), meanStandAge := 0]  ## NAs come from stands with 0 B and 0 age
vegDataForHVs[is.na(sdStandAge), sdStandAge := 0]  ## NAs come from stands with 0 B and 0 age or just one value of standAge

## calculate relative species B (across speciesCohorts)
vegDataForHVs[, standB := sum(B, na.rm = TRUE), by =  c("scenario", "rep", "year", "pixelIndex", "vegTypeCN")]
vegDataForHVs[, relB := sum(B) / standB, by = c("scenario", "rep", "year", "pixelIndex", "vegTypeCN", "speciesCode")]
vegDataForHVs[relB == "NaN", relB := 0]

## expand data
cols <- c("meanStandAge", "sdStandAge", "relB", "speciesCode", "scenario", "rep", "year", "pixelIndex", "vegTypeCN")   ## keep rep for wrapper.
vegDataForHVs <- unique(vegDataForHVs[, ..cols])
vegDataForHVs <- dcast.data.table(vegDataForHVs, as.formula("... ~ speciesCode"),
                                  value.var = "relB")
