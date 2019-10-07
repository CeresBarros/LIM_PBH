## clean workspace
rm(list=ls()); amc::.gc()

library(data.table)
library(ggplot2)

cohortDataFiles <- c(list.files("R/SpaDES/outputs/blogSep2019_noPM/", pattern = "cohortData", full.names = TRUE),
                     list.files("R/SpaDES/outputs/blogSep2019_PM/", pattern = "cohortData", full.names = TRUE))

allCohortData <- rbindlist(lapply(cohortDataFiles, FUN = function(ff) {
  cohortData <- readRDS(ff)
  yr <- as.integer(sub(".*year", "",  sub(".rds", "", ff)))
  scen <- if (grepl("_noPM", ff)) "noPM" else "PM"
  cohortData[, Year := yr]
  cohortData[, Scenario := scen]
  return(cohortData)
}))


ggplot(data = allCohortData,
       aes(x = Year, y = B, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  theme_classic() +
  facet_wrap(~ Scenario)

ggplot(data = allCohortData,
       aes(x = Year, y = B, fill = speciesCode)) +
  stat_summary(fun.y = mean, geom = "area", position = "fill") +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  theme_classic() +
  facet_wrap(~ Scenario)


ggplot(data = allCohortData,
       aes(x = Year, y = age, colour = speciesCode)) +
  stat_summary(fun.data = mean_sdl) +
  theme_classic() +
  facet_wrap(~ Scenario)



