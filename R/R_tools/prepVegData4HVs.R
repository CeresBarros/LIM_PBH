## VEGETATION DATA FOR HVs -----------------------

## Use the first fire year to identify the pixels we want to follow in time
## we follow the same pixels that were used to make fire attributes HVs, only now
## we select the start year and end years of the simulation
## the join shouldn't actually change anything because we already subset the pixels with veg
## in the montane belt (regardless of fire)

cacheExtra <- c(digest::digest(allPixelCohortDataMnt), digest::digest(summaryFireAttributes))
vegDataForHVs <- Cache(prepVegDataHVs,
                       allPixelCohortDataMnt = allPixelCohortDataMnt,
                       summaryFireAttributes = summaryFireAttributes,
                       useFirstLastYear = useFirstLastYear,
                       yearSubset = yearSubset,
                       yearSamples = yearSamples,
                       omitArgs = c("allPixelCohortDataMnt", "summaryFireAttributes"),
                       .cacheExtra = cacheExtra)

