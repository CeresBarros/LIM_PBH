## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)
library(raster)
library(quickPlot)
library(SpaDES)
simList_PM <- readRDS("R/SpaDES/outputs/blogSep2019_PM_oneFire/simList_fakeRstCurrentBurnblogSep2019_PM_oneFire.rds")
simList_noPM <- readRDS("R/SpaDES/outputs/blogSep2019_noPM_oneFire/simList_fakeRstCurrentBurnblogSep2019_noPM_oneFire.rds")

cohortDataFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire", pattern = "cohortData", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire", pattern = "cohortData", full.names = TRUE))

allCohortData <- rbindlist(fill = TRUE, l = lapply(cohortDataFiles, FUN = function(ff) {
  cohortData <- readRDS(ff)
  yr <- as.integer(sub(".*year", "",  sub(".rds", "", ff)))
  scen <- if (grepl("_noPM", ff)) "noPM" else "PM"
  cohortData[, Year := yr]
  cohortData[, Scenario := scen]
  return(cohortData)
}))

rstCurrentBurnFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "rstCurrentBurn", full.names = TRUE),
                         list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "rstCurrentBurn", full.names = TRUE))
pixelGroupMapFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "pixelGroupMap", full.names = TRUE),
                        list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "pixelGroupMap", full.names = TRUE))
vegTypeMapFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM_oneFire/", pattern = "vegTypeMap", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM_oneFire/", pattern = "vegTypeMap", full.names = TRUE))

rstCurrentBurnStk_noPM <- stack(lapply(grep("_noPM", rstCurrentBurnFiles, value = TRUE), readRDS))
rstCurrentBurnStk_PM <- stack(lapply(grep("_PM", rstCurrentBurnFiles, value = TRUE), readRDS))
pixelGroupMapStk_noPM <- stack(lapply(grep("_noPM", pixelGroupMapFiles, value = TRUE), readRDS))
pixelGroupMapStk_PM <- stack(lapply(grep("_PM", pixelGroupMapFiles, value = TRUE), readRDS))
vegTypeMapStk_noPM <- stack(lapply(grep("_noPM", vegTypeMapFiles, value = TRUE), readRDS))
vegTypeMapStk_PM <- stack(lapply(grep("_PM", vegTypeMapFiles, value = TRUE), readRDS))

names(rstCurrentBurnStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", rstCurrentBurnFiles, value = TRUE)))
names(rstCurrentBurnStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", rstCurrentBurnFiles, value = TRUE)))
names(pixelGroupMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", pixelGroupMapFiles, value = TRUE)))
names(pixelGroupMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", pixelGroupMapFiles, value = TRUE)))
names(vegTypeMapStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", vegTypeMapFiles, value = TRUE)))
names(vegTypeMapStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", vegTypeMapFiles, value = TRUE)))

## cheat

if (grepl("_oneFire", rstCurrentBurnFiles[1])) {
  ## compile all rasters except rstCurrentBurn
  pixelBurnCohortData <- rbind(
    rbindlist(lapply(names(vegTypeMapStk_noPM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(pixelGroupMapStk_noPM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
                 vegType = getValues(vegTypeMapStk_noPM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "noPM")
    })),
    rbindlist(lapply(names(vegTypeMapStk_PM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(pixelGroupMapStk_PM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
                 vegType = getValues(vegTypeMapStk_PM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "PM")
    }))
  )

  burnPix <- rbind(
    data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_noPM)),
               burnt = as.integer(!is.na(rstCurrentBurnStk_noPM[[1]][])),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_noPM))),
               Scenario = "noPM"),
    data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_PM)),
               burnt = as.integer(!is.na(rstCurrentBurnStk_noPM[[1]][])),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_PM))),
               Scenario = "PM")
  )

  ## checks
  if (length(setdiff(burnPix$pixelID, pixelBurnCohortData$pixelID)) |
      length(setdiff(pixelBurnCohortData$pixelID, burnPix$pixelID)))
    stop("pixel IDs differ between rstCurrentBurn and other output rasters")
  pixelBurnCohortData <- burnPix[, .(Scenario, pixelID, burnt)][pixelBurnCohortData, on = c("Scenario", "pixelID")]

  # pixelBurnCohortData <- na.omit(pixelBurnCohortData)

  ## checks
  testList <- split(pixelBurnCohortData[, .(burnt, Year)], by = "Year", keep.by = FALSE)
  testList <- lapply(testList, FUN = function(DT) DT[["burnt"]])

  test <- sapply(seq_along(testList), FUN =  function(n) {
    length(setdiff(testList[[n]], unlist(testList[-n])))
  })
  if (any(test > 0))
    stop("burnt pixel IDs differ between years")

} else {
  pixelBurnCohortData <- rbind(
    rbindlist(lapply(names(rstCurrentBurnStk_noPM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_noPM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_noPM[[x]]),
                 burnt = as.integer(!is.na(rstCurrentBurnStk_PM[[x]][])),
                 vegType = getValues(vegTypeMapStk_noPM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "noPM")
    })),
    rbindlist(lapply(names(rstCurrentBurnStk_PM), FUN = function(x) {
      data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_PM[[x]])),
                 pixelGroup = getValues(pixelGroupMapStk_PM[[x]]),
                 burnt = as.integer(!is.na(rstCurrentBurnStk_PM[[x]][])),
                 vegType = getValues(vegTypeMapStk_PM[[x]]),
                 Year = as.integer(sub("year", "", x)),
                 Scenario = "PM")
    }))
  )
  # pixelBurnCohortData <- na.omit(pixelBurnCohortData)
}

pixelBurnCohortData <- allCohortData[pixelBurnCohortData, nomatch = 0,
                                     on = c("Scenario", "Year", "pixelGroup")]

vegTypeTable <- as.data.table(levels(vegTypeMapStk_noPM[[1]]))
pixelBurnCohortData <- vegTypeTable[, .(ID, VALUE)][pixelBurnCohortData, on = "ID==vegType"]
setnames(pixelBurnCohortData, old = c("ID", "VALUE"),
         new = c("vegTypeID", "vegType"))

pixelBurnCohortData[, noPixels := length(pixelID), by = pixelGroup]

summaryBurnCohortData <- pixelBurnCohortData[, list(BiomassBySpecies = as.numeric(sum(B * noPixels, na.rm = TRUE)),
                                                    MortalityBySpecies = as.numeric(sum(mortality * noPixels, na.rm = TRUE)),
                                                    aNPPBySpecies = as.numeric(sum(aNPPAct * noPixels, na.rm = TRUE)),
                                                    AgeBySppWeighted = as.numeric(sum(age * B * noPixels, na.rm = TRUE) /
                                                                                    sum(B * noPixels, na.rm = TRUE))),
                                             by = .(Scenario, Year, burnt, speciesCode)]

speciesLabels <- c("Abie_sp" = "Fir", "Lari_lar" = "Larch",
                   "Pice_eng" = "En. spruce", "Pice_gla" = "Wh. spruce",
                   "Pice_mar" = "Bl. spruce", "Pinu_sp" = "Lo. pine",
                   "Popu_sp" = "Aspen", "Pseu_men" = "Douglas-fir")

plot1 <- ggplot(data = summaryBurnCohortData,
                aes(x = Year, y = log(BiomassBySpecies), colour = speciesCode)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "grey") +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2",
                      labels = speciesLabels) +
  labs(title = "Total landscape biomass", y = "log-Biomass (g/m^2)") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = summaryBurnCohortData,
                aes(x = Year, y = log(MortalityBySpecies), colour = speciesCode)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 5, linetype = "dashed", colour = "grey") +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  scale_colour_brewer(type = "qual", palette = "Dark2",
                      labels = speciesLabels) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot3 <- ggplot(data = pixelBurnCohortData,
                aes(x = Year, fill = vegType)) +
  geom_area(stat = "count", position = "fill") +
  scale_fill_brewer(type = "qual", palette = "Dark2",
                    labels = speciesLabels) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "No. pixels per vegetation type", y = "g/m^2") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

ggpubr::ggarrange(plot1, plot2, plot3, nrow = 3,
                  legend = "right", common.legend = TRUE)
ggsave(filename = "R/SpaDES/outputs/blogSep2019_noPM_PM_BMortVegType.tiff",
       width = 10, height = 15)


## make GIFs of vegetation maps

png(file="R/SpaDES/outputs/vegTypeMapStk_noPM%02d.png", width=200, height=200)
for (i in c(10:1, "G0!")){
  plot.new()
  text(.5, .5, i, cex = 6)
}
dev.off()

# convert the .png files to one .gif file using ImageMagick.
# The system() function executes the command as if it was done
# in the terminal. the -delay flag sets the time between showing
# the frames, i.e. the speed of the animation.
system("convert -delay 80 *.png example_1.gif")

# to not leave the directory with the single jpeg files
# I remove them.
file.remove(list.files(pattern=".png"))
Plot(vegTypeMapStk_PM)




## OTHER PLOTS
topoPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "BrBG"))
plot(simList_noPM$slopeRas, col = topoPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)

tempPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu"))
plot(simList_noPM$temperatureRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)

plot(simList_noPM$precipitationRas, col = tempPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)


sppPal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Blues"))
plot(simList_noPM$speciesLayers[[1]], col = sppPal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)

plot(simList_noPM$rawBiomassMap, axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)

agePal <- colorRampPalette(RColorBrewer::brewer.pal(n = 9, name = "Greens"))
plot(simList_noPM$standAgeMap, col = agePal(20), axes = FALSE,
     bg = "transparent", box = FALSE, legend = FALSE)


