## PLOT LABELS AND COLOURS --------------------------------

## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(summaryBurnCohortDataSpp$speciesCode), column = "EN_generic_full",
                                       df = preSimList$sppEquiv)
names(speciesLabels) <- unique(summaryBurnCohortDataSpp$speciesCode)

speciesColours <- levels(preSimList$vegTypeMap)[[1]]$colors
names(speciesColours) <- levels(preSimList$vegTypeMap)[[1]]$VALUE

## make vegType labels/colours
vegTypeLabels <- unique(levels(preSimList$vegTypeMap)[[1]]$VALUE)
vegTypeLabels <- LandR::equivalentName(value = vegTypeLabels, column = "EN_generic_full",
                                       df = preSimList$sppEquiv)
vegTypeLabels[length(vegTypeLabels) + 1] <- "No veg."
names(vegTypeLabels) <- c(levels(preSimList$vegTypeMap)[[1]]$ID, "0")

vegTypeColours <- as.character(levels(preSimList$vegTypeMap)[[1]]$colors)
vegTypeColours[length(vegTypeColours) + 1] <- "grey40"
names(vegTypeColours) <- c(levels(preSimList$vegTypeMap)[[1]]$ID, "0")


vegTypeCNLabels <- unique(as.character(allHVData$vegType))
names(vegTypeCNLabels) <- vegTypeCNLabels
vegTypeCNLabels <- sub("PIEN", "Spruce", vegTypeCNLabels)
vegTypeCNLabels <- sub("MMC", "Moist conif.", vegTypeCNLabels)
vegTypeCNLabels <- sub("mixedwood", "Mixed", vegTypeCNLabels)
vegTypeCNLabels <- sub("broadleaf", "Broadleaf", vegTypeCNLabels)
vegTypeCNLabels <- sub("PICO", "Pine", vegTypeCNLabels)
vegTypeCNLabels <- sub("DMCPSME", "Dry conif./Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("^PSME$", "Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("dryPSME", "Dry Douglas-fir", vegTypeCNLabels)

## reorder
vegTypeCNLabels <- vegTypeCNLabels[c(grep("No veg.|landscape", vegTypeCNLabels, invert = TRUE),
                                     grep("No veg.", vegTypeCNLabels),
                                     grep("landscape", vegTypeCNLabels))]

## landscape gets a different colour
vegTypeCNColours <- vegTypeCNLabels
vegTypeCNColours[1:length(vegTypeCNColours)-1] <- RColorBrewer::brewer.pal(length(vegTypeCNColours)-1, name = "Set1")
vegTypeCNColours["landscape"] <- "darkgreen"
