## FIRE DATA SUMMARY FOR HVs -----------------------
## summarize data first - as in Steel et al 2021, fire properties are summarized across time, but by pixel
## (and by scenario/rep)
## only look at pixels with vegetation dynamics so that we can compare with biodiv. HVs
cacheExtra <- c(digest::digest(allPixelBurnData), digest::digest(allPixelCohortDataMnt))
summaryFireAttributes <- Cache(summariseFireRegAttrs,
                               allPixelBurnData = allPixelBurnData,
                               allPixelCohortDataMnt = allPixelCohortDataMnt,
                               yearSubset = yearSubset,
                               omitArgs = c("allPixelBurnData", "allPixelCohortDataMnt"),
                               .cacheExtra = cacheExtra)
