library(data.table)
library(raster)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(LandR)

files <- list.files(file.path("R/SpaDES/outputs/mar2022Runs/noPM"), "cohortData", full.names = TRUE)
files <- grep(paste("year", seq(2011, 2111, 5), sep = "", collapse = "|"), files, value = TRUE)

files2 <- list.files(file.path("R/SpaDES/outputs/mar2022Runs/noPM"), "pixelGroupMap", full.names = TRUE)
files2 <- grep(paste("year", seq(2011, 2111, 5), sep = "", collapse = "|"), files2, value = TRUE)

allCohortDataMarch <- mapply(cdfile = files, pgfile = files2, FUN = function(cdfile, pgfile) {
  if (as.integer(sub("\\.rds", "", sub(".*year", "", cdfile))) ==
      as.integer(sub("\\.rds", "", sub(".*year", "", pgfile)))) {
    dat <- readRDS(cdfile)
    map <- readRDS(pgfile)
    dat$year <- as.integer(sub("\\.rds", "", sub(".*year", "", cdfile)))
    dat <- addPixels2CohortData(cohortData = dat, pixelGroupMap = map)
    return(dat)
  } else stop("file years don't match")
}, SIMPLIFY = FALSE)
allCohortDataMarch <- rbindlist(allCohortDataMarch, use.names = TRUE)

files <- list.files(file.path("R/SpaDES/outputs/mar2022Runs/noPM"), "severityData", full.names = TRUE)
files <- grep(paste("year", seq(2011, 2111, 1), sep = "", collapse = "|"), files, value = TRUE)

allSevDataMarch <- lapply(files, function(ff) {
  year <- as.integer(sub("\\.rds", "", sub(".*year", "", ff)))
  dat <- readRDS(ff)
  dat$year <- year
  dat
})
allSevDataMarch <- rbindlist(allSevDataMarch, use.names = TRUE)
allSevDataMarch[, year := year - 1] ## cohortData is saved after fire so the lost biomass reflects last years total.

allCohortDataMarch[, totalB := sum(B), by = .(pixelIndex, year)]
allCohortDataMarch[, totalMort := sum(mortality), by = .(pixelIndex, year)]
allCohortDataMarch[, totalANPP := sum(aNPPAct), by = .(pixelIndex, year)]
testMar <- allSevDataMarch[allCohortDataMarch, on = .(pixelIndex, year)]
if (nrow(testMar[!is.na(severityB)][severityB != totalB])) {
  stop("this years burnt B must be equal to last year's total")
}

files <- unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", list = TRUE)$Name
files <- list(files = grep("noPM/noPM_rep01/cohortData", files, value = TRUE),
              files2 = grep("noPM/noPM_rep01/pixelGroupMap", files, value = TRUE))
files <- lapply(files, function(files) unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", files = files, exdir = "F:/LandscapesInMotion/R/SpaDES/outputs"))
files2 <- files$files2
files <- files$files

allCohortDataJune <- mapply(cdfile = files, pgfile = files2, FUN = function(cdfile, pgfile) {
  if (as.integer(sub("\\.rds", "", sub(".*year", "", cdfile))) ==
      as.integer(sub("\\.rds", "", sub(".*year", "", pgfile)))) {
    dat <- readRDS(cdfile)
    map <- readRDS(pgfile)
    dat$year <- as.integer(sub("\\.rds", "", sub(".*year", "", cdfile)))
    dat <- addPixels2CohortData(cohortData = dat, pixelGroupMap = map)
    return(dat)
  } else stop("file years don't match")
}, SIMPLIFY = FALSE)
allCohortDataJune <- rbindlist(allCohortDataJune, use.names = TRUE)

files <- unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", list = TRUE)$Name
files <- grep("noPM/noPM_rep01/severityData", files, value = TRUE)
unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", files = files, exdir = "F:/LandscapesInMotion/R/SpaDES/outputs")
files <- file.path("F:/LandscapesInMotion/R/SpaDES/outputs", files)

allSevDataJune <- lapply(files, function(ff) {
  year <- as.integer(sub("\\.rds", "", sub(".*year", "", ff)))
  dat <- readRDS(ff)
  dat$year <- year
  dat
})
allSevDataJune <- rbindlist(allSevDataJune, use.names = TRUE)
allSevDataJune[, year := year - 1] ## cohortData is saved after fire so for the lost biomass reflects last years total.

allCohortDataJune[, totalB := sum(B), by = .(pixelIndex, year)]
allCohortDataJune[, totalMort := sum(mortality), by = .(pixelIndex, year)]
allCohortDataJune[, totalANPP := sum(aNPPAct), by = .(pixelIndex, year)]
testJun <- allSevDataJune[allCohortDataJune, on = .(pixelIndex, year)]
if (nrow(testJun[!is.na(severityB)][severityB != totalB])) {
  stop("this years burnt B must be equal to last year's total")
}

plotData1 <- melt(testJun, measure.vars = c("totalB", "totalMort", "totalANPP"),
                  variable.name = "variable", value.name = "B")
plotData2 <- melt(testMar, measure.vars = c("totalB", "totalMort", "totalANPP"),
                  variable.name = "variable", value.name = "B")
limits <- c(min(plotData1[, sum(B.1, na.rm = TRUE), by = .(year, variable)]$V1, plotData2[, sum(B.1, na.rm = TRUE), by = .(year, variable)]$V1),
            max(plotData1[, sum(B.1, na.rm = TRUE), by = .(year, variable)]$V1, plotData2[, sum(B.1, na.rm = TRUE), by = .(year, variable)]$V1))

plot1 <- ggplot(plotData1, aes(x = year, y = B.1, colour = variable)) +
  stat_summary(fun = "sum", geom = "line", size = 1) +
  stat_summary_bin(data = allSevDataJune, mapping = aes(x = year, y = severityB, col = "severityB"),
                   fun = "sum", geom = "line", size = 1, breaks = seq(2011, 2111, 5)) +
  geom_vline(xintercept = min(allSevDataJune$year), linetype = "dashed", colour = "black", size = 0.7) +
  geom_vline(xintercept = seq(2011, 2111, 10), linetype = "dotted", colour = "grey", size = 0.3) +
  scale_y_continuous(limits = limits) +
  scale_x_continuous(breaks = unique(allSevDataJune$year)) +
  theme_classic() +
  labs(title = "old results")

plot2 <- ggplot(plotData2, aes(x = year, y = B.1, colour = variable)) +
  stat_summary(fun = "sum", geom = "line", size = 1) +
  stat_summary_bin(data = allSevDataMarch, mapping = aes(x = year, y = severityB, col = "severityB"),
                   fun = "sum", geom = "line", size = 1, breaks = seq(2011, 2111, 5)) +
  # stat_summary(data = allSevDataMarch, mapping = aes(x = year, y = severityB, col = "severityB"),
  # fun = "sum", geom = "line", size = 1) +
  geom_vline(xintercept = min(allSevDataMarch$year), linetype = "dashed", colour = "black", size = 0.7) +
  geom_vline(xintercept = seq(2011, 2111, 10), linetype = "dotted", colour = "grey", size = 0.3) +
  # scale_y_continuous(limits = limits) +
  scale_x_continuous(breaks = unique(allSevDataMarch$year)) +
  theme_classic() +
  labs(title = "new results")

ggarrange(plot1, plot2, nrow = 1)

# simInitJune <- grep("LIM_simInit_noPM", unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", list = TRUE), value = TRUE)
# simInitJune <- unzip("F:/LandscapesInMotion/R/SpaDES/outputs/jun2021Runs.zip", files = simInitJune)
# SpaDES.core::loadSimList("F:/LandscapesInMotion/R/SpaDES/R/SpaDES/outputs/mar2022Runs/noPM/LIM_simInit_noPM.qs")
# simInitMarch <- SpaDES.core::loadSimList("F:/LandscapesInMotion/R/SpaDES/R/SpaDES/outputs/mar2022Runs/noPM/LIM_simInit_noPM.qs")
# firstYearB <- allCohortDataMarch[year == 2011]

# files <- list.files("R/SpaDES/outputs/mar2022Runs/noPM/", "rstCurrentFires", full.names = TRUE)
# files <- list.files("R/SpaDES/outputs/mar2022Runs/noPM/noPM_rep1/", "rstCurrentFires", full.names = TRUE)
# files <- list.files("R/SpaDES/outputs/mar2022Runs/PM/", "rstCurrentFires", full.names = TRUE)
files <- list.files("R/SpaDES/outputs/mar2022Runs/PM/PM_rep1/", "rstCurrentFires", full.names = TRUE)
# files <- grep(paste("year", seq(2011, 2111, 1), sep = "", collapse = "|"), files, value = TRUE)

rstCurrentFires <- stack(lapply(files, readRDS))
names(rstCurrentFires) <- sub("\\.rds", "", sub(".*year", "", files))

fireSizes <- lapply(unstack(rstCurrentFires), function(ras) {
  tab <- table(ras[])
  if (dim(tab)) {
    as.data.table(tab)
  } else NULL
})
fireSizes <- rbindlist(fireSizes)
setnames(fireSizes, new = c("fireID", "noPixels"))

fireRasterObs <- reproducible::Cache(LandR::prepInputsFireYear,
                                     url = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip",
                                     destinationPath = simPaths$inputPath,
                                     rasterToMatch = LIM_simInitList$PM$rasterToMatch,
                                     maskWithRTM = TRUE,
                                     method = "ngb",
                                     datatype = "INT2U",
                                     filename2 = NULL,
                                     fireField = "CFS_REF_ID",
                                     earliestYear = 1,
                                     fun = "sf::st_read",
                                     userTags = c("FavierFireSpread", "function:.inputObjects"),
                                     omitArgs = c("destinationPath", "targetFile", "userTags"))
fireSizesObs <- as.data.table(table(fireRasterObs[]))
setnames(fireSizesObs, new = c("fireID", "noPixels"))

allFireSizes <- rbindlist(list("simulated" = fireSizes, "observed" = fireSizesObs), idcol = TRUE)

ggplot2::ggplot(allFireSizes, ggplot2::aes(x = noPixels, fill = .id)) +
  ggplot2::geom_density(alpha = 0.5)

