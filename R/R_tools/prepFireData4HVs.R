## FIRE DATA SUMMARY FOR HVs -----------------------
## summarize data first - as in Steel et al 2021, fire properties are summarized across time, but by pixel
## (and by scenario/rep)
## only look at pixels with vegetation dynamics so that we can compare with biodiv. HVs
summaryFireAttributes <- allPixelBurnData[!is.na(pixelGroup), list(meanFreq = mean(fireFreq),
                                                                   meanSev = mean(severity),
                                                                   meanSevB = mean(severityB),
                                                                   meanPatchS = mean(patchSize)),
                                          by = .(scenario, rep, pixelIndex)]
## add vegType per pixel at the start of the simulation
## not the first fire year, as these vegTypes have already been altered by fire
## and add pixels that had no fires
cols <- c("pixelIndex", "vegTypeCN", "scenario", "rep")
summaryFireAttributes <- summaryFireAttributes[unique(allPixelCohortDataMnt[year == max(yearSubset), ..cols]),
                                               on = c("scenario", "rep", "pixelIndex")]
## checks
test1 <- any(is.na(summaryFireAttributes$vegTypeCN))
if (test1) {
  stop("NA vegTypeCNs where pixelGroup - i.e. vegetation - exists")
}

test2 <- allPixelBurnData[!is.na(pixelGroup)][summaryFireAttributes[is.na(meanFreq), .(scenario, rep, pixelIndex)],
                                              on = .(scenario, rep, pixelIndex), nomatch = 0]
if (dim(test2)[1]) {
  stop("pixels that had fire and veg data in allPixelBurnData were ",
       "accidentally dropped when adding vegTypeCN")
}

## set fire attributes to 0 in pixels that had no fires
cols <- c("meanFreq", "meanSev", "meanSevB", "meanPatchS")
summaryFireAttributes[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]
amc::.gc()
