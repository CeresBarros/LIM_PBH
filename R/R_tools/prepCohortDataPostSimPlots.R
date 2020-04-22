## assuming you have a cohortData table, pixelGroupMap and vegTypeMap
## all from the save year and saved at the same point in a timestep

## Step 1. Add pixels to cohortData
pixelCohortData <- addPixels2CohortData(cohortData, pixelGroupMap)

## Step 2. Make a table of vegType data
vegTypeData <- data.table(pixelIndex = seq_len(ncell(vegTypeMap)),
                          pixelGroup = getValues(pixelGroupMap),
                          vegType = vegTypeMap)

vegTypeData <- vegTypeData[!is.na(pixelGroup)]

## Step 3. Join
## note that vegTypeData and pixelBurnData have more pixels because PGs of 0 are there, but not on cohortData
## so join by keeping all pixels
pixelCohortData <- merge(pixelCohortData, vegTypeData,
                         by = c("pixelIndex", "pixelGroup"), all = TRUE)

## Step 4 ADD MISSING SPECIES IN PIXEL COMBINATION
## cohortData doens't track absent cohorts, so they need to ba added back
## for now pixels from the 0s pixelGroup  have one entry with NAs for speciesCode
## they will be ignored for now and removed later, after adding one species entry for each of these pixels.
## for reporting consistency add to show losses in B
combinations <- unique(pixelCohortData[, .(pixelIndex, pixelGroup)])
spp <- as.character(na.omit(unique(pixelCohortData$speciesCode)))
combinations <- lapply(spp, FUN = function(x) {
  data.table(combinations,
             speciesCode = x)
}) %>%
  rbindlist(., use.names = TRUE)

## join while keeping all combos, NA species will now disappear.
pixelCohortData <- pixelCohortData[combinations, on = .(pixelIndex, pixelGroup, speciesCode)]

## add ecoregion group where it's missing
pixelCohortData[, ecoregionGroup := unique(na.omit(ecoregionGroup)),
                by = .(pixelGroup)]
## add vegType where it's missing, but it's a pixel with some veg
pixelCohortData[, vegType := max(vegType, na.rm = TRUE),
                by = .(pixelGroup)]

## replace NAs of cohortData by 0s and add missing spp.
cols <- c("age", "B", "mortality", "aNPPAct", "vegType")
replaceNAs <- function(x, val = 0) {
  x[is.na(x)] <- val
  x
}

pixelCohortData[, (cols) := lapply(.SD, replaceNAs), .SDcols = cols]

## Step 5 SUMMARY ACROSS LANDSCAPE - maybe you don't need this?
## remember that biomass is multiplied by 100 in *boreal*, this will revert the units to tonnes/ha
# summaryCohortData <- pixelCohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
#                                             MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
#                                             aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
#                                             AgeBySppWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
#                                                                             sum((B/100), na.rm = TRUE))),
#                                      by = .(speciesCode)]
# summaryCohortData[is.nan(AgeBySppWeighted), AgeBySppWeighted := 0]

## Step 6. make labels
## make species labels/colours
speciesLabels <- LandR::equivalentName(value = unique(pixelCohortData$speciesCode), column = "EN_generic_full",
                                       df = sppEquiv)
names(speciesLabels) <- unique(pixelCohortData$speciesCode)

speciesColours <- levels(vegTypeMap[[1]])[[1]]$colors
names(speciesColours) <- levels(vegTypeMap[[1]])[[1]]$VALUE

## make vegType labels/colours
vegTypeLabels <- as.character(levels(vegTypeMap[[1]])[[1]]$VALUE)
vegTypeLabels <- LandR::equivalentName(value = vegTypeLabels, column = "EN_generic_full",
                                       df = sppEquiv)
vegTypeLabels[length(vegTypeLabels) + 1] <- "No veg."
names(vegTypeLabels) <- c(levels(vegTypeMap[[1]])[[1]]$ID, "0")

vegTypeColours <- as.character(levels(vegTypeMap[[1]])[[1]]$colors)
vegTypeColours[length(vegTypeColours) + 1] <- "grey40"
names(vegTypeColours) <- c(levels(vegTypeMap[[1]])[[1]]$ID, "0")

## Step 7. Plots
plotData <- pixelCohortData[, list(noPixelsVeg = length(unique(pixelIndex))),
                               by = .(scenario, year, firePresAbs, vegType)]
plot1 <- ggplot(data = plotData,
                aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "stack") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "no. pixels") +


plot2 <- ggplot(data = plotData,
                  aes(x = year, y = noPixelsVeg, fill = as.factor(vegType))) +
  geom_area(stat = "identity", position = "fill") +
  theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = vegTypeColours, labels = vegTypeLabels) +
  labs(title = "Dominant species", y = "prop. pixels")


