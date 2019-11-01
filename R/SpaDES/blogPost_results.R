## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)
library(raster)
library(quickPlot)
simList_PM <- readRDS("R/SpaDES/outputs/blogSep2019_PM/simList_fakeRstCurrentBurnblogSep2019_PM.rds")
simList_noPM <- readRDS("R/SpaDES/outputs/blogSep2019_noPM/simList_fakeRstCurrentBurnblogSep2019_noPM.rds")   ## NEEDS TO BE REDONE

cohortDataFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM/", pattern = "cohortData", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM/", pattern = "cohortData", full.names = TRUE))

allCohortData <- rbindlist(fill = TRUE, l = lapply(cohortDataFiles, FUN = function(ff) {
  cohortData <- readRDS(ff)
  yr <- as.integer(sub(".*year", "",  sub(".rds", "", ff)))
  scen <- if (grepl("_noPM", ff)) "noPM" else "PM"
  cohortData[, Year := yr]
  cohortData[, Scenario := scen]
  return(cohortData)
}))

rstCurrentBurnFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM/", pattern = "rstCurrentBurn", full.names = TRUE),
                         list.files("R/SpaDES/outputs/blogSep2019_PM/", pattern = "rstCurrentBurn", full.names = TRUE))

rstCurrentBurnStk_noPM <- stack(lapply(grep("_noPM", rstCurrentBurnFiles, value = TRUE), readRDS))
rstCurrentBurnStk_PM <- stack(lapply(grep("_PM", rstCurrentBurnFiles, value = TRUE), readRDS))

names(rstCurrentBurnStk_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", rstCurrentBurnFiles, value = TRUE)))
names(rstCurrentBurnStk_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", rstCurrentBurnFiles, value = TRUE)))

plot(rstCurrentBurnStk_noPM)
plot(rstCurrentBurnStk_PM)

pixelGroupMapFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM/", pattern = "pixelGroupMap", full.names = TRUE),
                        list.files("R/SpaDES/outputs/blogSep2019_PM/", pattern = "pixelGroupMap", full.names = TRUE))

pixelGroupMap_noPM <- stack(lapply(grep("_noPM", pixelGroupMapFiles, value = TRUE), readRDS))
pixelGroupMap_PM <- stack(lapply(grep("_PM", pixelGroupMapFiles, value = TRUE), readRDS))

names(pixelGroupMap_noPM) <- sub(".*year", "year", sub("\\.rds", "", grep("_noPM", pixelGroupMapFiles, value = TRUE)))
names(pixelGroupMap_PM) <- sub(".*year", "year", sub("\\.rds", "", grep("_PM", pixelGroupMapFiles, value = TRUE)))

plot(pixelGroupMap_noPM)
plot(pixelGroupMap_PM)

burnPixelTable <- rbind(
  rbindlist(lapply(names(rstCurrentBurnStk_noPM),
                   FUN = function(x) {
                     data.table(pixID = seq_len(ncell(rstCurrentBurnStk_noPM[[x]])),
                                pixelGroup = getValues(pixelGroupMap_noPM[[x]]),
                                burnt = as.integer(!is.na(getValues(rstCurrentBurnStk_noPM[[x]]))),
                                Year = as.integer(sub("year", "", x)),
                                Scenario = "noPM")
                   })),
  rbindlist(lapply(names(rstCurrentBurnStk_PM),
                   FUN = function(x) {
                     data.table(pixID = seq_len(ncell(rstCurrentBurnStk_PM[[x]])),
                                pixelGroup = getValues(pixelGroupMap_PM[[x]]),
                                burnt = as.integer(!is.na(getValues(rstCurrentBurnStk_PM[[x]]))),
                                Year = as.integer(sub("year", "", x)),
                                Scenario = "PM")
                   }))
)

burnPixelTable <- na.omit(burnPixelTable)

setkey(burnPixelTable, Scenario, Year, pixelGroup)
setkey(allCohortData, Scenario, Year, pixelGroup)

## add cohort data to burnt pixels
burnPixelCohortData <- allCohortData[burnPixelTable]

plot1 <- ggplot(data = burnPixelCohortData,
                aes(x = Year, y = B, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  stat_summary(fun.y = mean, geom = "line") +
  scale_colour_brewer(type = "qual", palette = "Set3") +
  theme_classic() +
  theme(strip.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "No fire", "1" = "Fire"),)) +
  labs(y = "Biomass g/m^2")


plot2 <- ggplot(data = burnPixelCohortData,
                aes(x = Year, y = mortality, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  stat_summary(fun.y = mean, geom = "line") +
  scale_colour_brewer(type = "qual", palette = "Set3") +
  theme_classic() +
  theme(strip.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.text = element_text(size = 16)) +
  facet_grid(burnt ~ Scenario,
             labeller = labeller(burnt = c("0" = "No fire", "1" = "Fire"))) +
  labs(y = "Mortality (biomass loss, g/m^2)")

ggpubr::ggarrange(plot1, plot2, common.legend = TRUE, legend = "bottom")
