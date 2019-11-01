## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)
library(raster)
library(quickPlot)
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

names(rstCurrentBurnStk_noPM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_noPM", rstCurrentBurnFiles, value = TRUE)))
names(rstCurrentBurnStk_PM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_PM", rstCurrentBurnFiles, value = TRUE)))
names(pixelGroupMapStk_noPM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_noPM", pixelGroupMapFiles, value = TRUE)))
names(pixelGroupMapStk_PM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_PM", pixelGroupMapFiles, value = TRUE)))
names(vegTypeMapStk_noPM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_noPM", vegTypeMapFiles, value = TRUE)))
names(vegTypeMapStk_PM) <- sub(".*year(0)*", "year", sub("\\.rds", "", grep("_PM", vegTypeMapFiles, value = TRUE)))

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
               burnt = getValues(rstCurrentBurnStk_noPM[[1]]),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_noPM))),
               Scenario = "noPM"),
    data.table(pixelID = seq_len(ncell(rstCurrentBurnStk_PM)),
               burnt = getValues(rstCurrentBurnStk_PM[[1]]),
               Year = as.integer(sub("year", "", names(rstCurrentBurnStk_PM))),
               Scenario = "PM")
    )

  ## checks
  if (length(setdiff(burnPix$pixelID, pixelBurnCohortData$pixelID)) |
      length(setdiff(pixelBurnCohortData$pixelID, burnPix$pixelID)))
    stop("pixel IDs differ between rstCurrentBurn and other output rasters")
  pixelBurnCohortData <- burnPix[, .(Scenario, pixelID, burnt)][pixelBurnCohortData, on = c("Scenario", "pixelID")]

  pixelBurnCohortData <- na.omit(pixelBurnCohortData)

  ## checks
  outer(X = lapply(split(pixelBurnCohortData, f = "Year"),
                   FUN = function(DT) DT[burnt == 1, pixelID]),
        Y = lapply(split(pixelBurnCohortData, f = "Year"),
                   FUN = function(DT) DT[burnt == 1, pixelID]),
        setdiff)

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
}

pixelBurnCohortData <- na.omit(pixelBurnCohortData)
pixelBurnCohortData <- allCohortData[pixelBurnCohortData, on = c("Scenario", "Year", "pixelGroup")]

vegTypeTable <- as.data.table(levels(vegTypeMapStk_noPM[[1]]))
pixelBurnCohortData <- vegTypeTable[, .(ID, VALUE)][pixelBurnCohortData, on = "ID==vegType"]
setnames(pixelBurnCohortData, old = c("ID", "VALUE"),
         new = c("vegTypeID", "vegType"))

plot1 <- ggplot(data = pixelBurnCohortData,
       aes(x = Year, y = B, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "Total across-landscape biomass", y = "g/m^2") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot2 <- ggplot(data = pixelBurnCohortData,
       aes(x = Year, y = mortality, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  stat_summary(fun.y = mean, geom = "line", size = 1) +
  theme_classic() +
  theme(text = element_text(size = 16),
        legend.title = element_blank()) +
  labs(title = "Total across-landscape mortality", y = "g/m^2") +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "no fire", "1" = "fire")))

plot3 <- ggplot(data = pixelBurnCohortData,
       aes(x = Year, colour = vegType)) +
  geom_bar(stat = "count", position = "dodge") +
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
