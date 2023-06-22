## FIRE DATA SUMMARY FOR HVs -----------------------
## summarize data first - as in Steel et al 2021, fire properties are summarized across time, but by pixel
## (and by scenario/rep)
## only look at pixels with vegetation dynamics so that we can compare with biodiv. HVs

summaryFireAttributes <- allPixelBurnData[pixelIndex %in% allPixelCohortDataMnt$pixelIndex]
summaryFireAttributes <- summaryFireAttributes[, list(meanFreq = unique(fireFreq),   ## note that fireFreq is already an average of mean fire intervals
                                                      meanSev = mean(severity),
                                                      meanSevB = mean(severityB),
                                                      meanPatchS = mean(patchSizeLogHa)),
                                               by = .(scenario, rep, pixelIndex)]

## add vegType per pixel at the end of the simulation
## and add pixels that had no fires
cols <- c("pixelIndex", "vegTypeCN", "scenario", "rep")
summaryFireAttributes <- summaryFireAttributes[unique(allPixelCohortDataMnt[year == max(yearSubset), ..cols]),
                                               on = c("scenario", "rep", "pixelIndex")]
## checks
if (getOption("LandR.assertions", TRUE)) {
  test1 <- any(is.na(summaryFireAttributes$vegTypeCN))
  if (test1) {
    stop("NA vegTypeCNs where pixelGroup - i.e. vegetation - exists")
  }

  test2 <- allPixelBurnData[pixelIndex %in% allPixelCohortDataMnt$pixelIndex][summaryFireAttributes[is.na(meanFreq), .(scenario, rep, pixelIndex)],
                                                                              on = .(scenario, rep, pixelIndex), nomatch = 0]
  if (nrows(test2)) {
    stop("pixels that had fire and veg data in allPixelBurnData were ",
         "accidentally dropped when adding vegTypeCN")
  }

  test3 <- any(is.na(summaryFireAttributes$meanFreq))
  test4 <- any(is.na(summaryFireAttributes$meanSev))
  test5 <- any(is.na(summaryFireAttributes$meanSevB))
  test6 <- any(is.na(summaryFireAttributes$meanPatchS))

  if (any(test3, test4, test5, test6)) {
    stop("Found NAs in fire properties")
  }

  rm(test1, test2, test3, test4, test5, test6)
}
gc(reset = TRUE)
