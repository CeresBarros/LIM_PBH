## PLOT LABELS AND COLOURS --------------------------------

## make species labels/colours
if (exists("summaryBurnCohortDataSpp")) {
  allSpeciesCodes <- unique(summaryBurnCohortDataSpp$speciesCode)
} else {
  if (exists("allPixelCohortDataMnt")) {
    allSpeciesCodes <- unique(allPixelCohortDataMnt$speciesCode)
  }
}

if (exists("allSpeciesCodes")) {
  sppEquiv <- preSimList$sppEquiv
  sppEquiv[, Latin_short := sub("(^[[:alpha:]]*) ([[:alpha:]]{1,3}).*", "\\1 \\2.", Latin_full)]
  speciesLabels <- LandR::equivalentName(value = allSpeciesCodes, column = "Latin_short",
                                         df = sppEquiv)
  names(speciesLabels) <- allSpeciesCodes
  speciesLabels[names(speciesLabels) == "Popu_sp"] <- "Populus/Betula spp"
  speciesLabels[names(speciesLabels) == "Pinu_sp"] <- "Pinus spp"
  speciesLabels[names(speciesLabels) == "Abie_sp"] <- "Abies spp"

  speciesColours <- levels(preSimList$vegTypeMap)[[1]]$colors
  names(speciesColours) <- levels(preSimList$vegTypeMap)[[1]]$VALUE
}

## make vegType labels/colours
vegTypeLabels <- unique(levels(preSimList$vegTypeMap)[[1]]$VALUE)
vegTypeLabels <- LandR::equivalentName(value = vegTypeLabels, column = "EN_generic_full",
                                       df = preSimList$sppEquiv)
vegTypeLabels[length(vegTypeLabels) + 1] <- "No veg."
names(vegTypeLabels) <- c(levels(preSimList$vegTypeMap)[[1]]$ID, "0")

vegTypeLabels[vegTypeLabels == "Deciduous"] <- "Broadleaf"
vegTypeLabels <- sort(vegTypeLabels)

vegTypeColours <- as.character(levels(preSimList$vegTypeMap)[[1]]$colors)
vegTypeColours[length(vegTypeColours) + 1] <- "grey40"
names(vegTypeColours) <- c(levels(preSimList$vegTypeMap)[[1]]$ID, "0")

if (exists("allPixelCohortDataMnt")) {
  vegTypeCNLabels <- unique(as.character(allPixelCohortDataMnt$vegTypeCN))
} else {
  vegTypeCNLabels <- unique(as.character(allHVData$vegType))
}

vegTypeCNLabels <- sort(vegTypeCNLabels)

names(vegTypeCNLabels) <- vegTypeCNLabels
vegTypeCNLabels <- sub("PIEN", "Spruce", vegTypeCNLabels)
vegTypeCNLabels <- sub("MMC", "Moist conif.", vegTypeCNLabels)
vegTypeCNLabels <- sub("mixedwood", "Mixed", vegTypeCNLabels)
vegTypeCNLabels <- sub("broadleaf", "Broadleaf", vegTypeCNLabels)
vegTypeCNLabels <- sub("PICO", "Pine", vegTypeCNLabels)
vegTypeCNLabels <- sub("DMCPSME", "Dry conif.", vegTypeCNLabels)
vegTypeCNLabels <- sub("^PSME$", "Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("dryPSME", "Dry Douglas-fir", vegTypeCNLabels)
vegTypeCNLabels <- sub("No veg.", "No forest", vegTypeCNLabels)
vegTypeCNLabels["landscape"] <- "landscape"

## reorder
vegTypeCNLabels <- vegTypeCNLabels[c(grep("No forest|landscape", vegTypeCNLabels, invert = TRUE),
                                     grep("No forest", vegTypeCNLabels),
                                     grep("landscape", vegTypeCNLabels))]

## landscape gets a different colour
vegTypeCNColours <- RColorBrewer::brewer.pal(length(vegTypeCNLabels) - 2, name = "Dark2")
## replace last colour (grey) by something else otherwise it's the same as no veg
if (length(vegTypeCNColours) == 8) {
  vegTypeCNColours[8] <- RColorBrewer::brewer.pal(3, name = "Paired")[2]
}

names(vegTypeCNColours) <- names(vegTypeCNLabels)[1:(length(vegTypeCNLabels) - 2)]
vegTypeCNColours["landscape"] <- "darkgreen"
vegTypeCNColours["No veg."] <- "grey40"


## labels for scenario:
scenLabels <- c("PM" = expression(M[MS]), "noPM" = expression(M[SR]),
                "HV_PM" = expression(M[MS]), "HV_noPM" = expression(M[SR]))
