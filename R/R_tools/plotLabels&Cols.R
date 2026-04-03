## PLOT LABELS AND COLOURS --------------------------------

if (!exists("preSimList")) {
  preSimList <- loadSimList(file.path(simPaths$outputPath, "LIM_preSimulation.qs"))
}

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
  speciesLabels <- sort(speciesLabels)
  speciesLabels["landscape"] <- "landscape"


  speciesColours <- levels(preSimList$vegTypeMap)[[1]]$colors
  names(speciesColours) <- levels(preSimList$vegTypeMap)[[1]]$VALUE
  speciesColours["landscape"] <- "#006400"
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
vegTypeColours[length(vegTypeColours) + 1] <- "#666666"
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
vegTypeCNColours["landscape"] <- "#006400"
vegTypeCNColours["No veg."] <- "#666666"

## colours for scenario, with "observed"
scenColours <- c("noPM" = scales::hue_pal()(2)[1], "PM" = scales::hue_pal()(2)[2],
                 "observed" = "grey50")

## labels for scenario:
scenLabels <- c("PM" = expression(M[MS]), "noPM" = expression(M[SR]), "observed" = "observed",
                "HV_PM" = expression(M[MS]), "HV_noPM" = expression(M[SR]))

## line type for scenario
scenLinetype <- c("PM" = 1, "noPM" = 2,
                  "HV_PM" = 1, "HV_noPM" = 1)

## fire labels
fireLabels <- c("0" = "no fire", "1" = "fire")
