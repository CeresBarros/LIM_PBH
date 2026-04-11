## GET CAMERON'S AGE DATA AND STAND VEG TYPES -----------------------
source("data/CameronsAgeData/Stand age function_for Ceres_Adapted.R")
ageDataCN <- fread("data/CameronsAgeData/treelist_outputs_for Ceres v2.csv")   ## v2 has stand density and basal area
patchVegTypeCN <- fread("data/CameronsAgeData/patch outputs_for Ceres_Oct2021.csv")

## remove records after 1940 (fire exclusion)
ageDataCN <- ageDataCN[Establishment.date <= 1940]

## remove funky record
ageDataCN <- ageDataCN[Reconstructed.age != 2018]

## calculate basal-area-weighted stand-ages for FT-level comparisons
cols <- c("Stand.age.1", "Stand.age.2", "Stand.age.3", "Species.age.3")
ageDataCN[, (cols) := stand.age.func(.SD), by = Patch.ID]

## replace NAs with stand age.
cols <- c("Stand.age.1", "Stand.age.2", "Stand.age.3")
ageDataCN[, (cols) := lapply(.SD, function(x) unique(x[!is.na(x)])), by = Patch.ID, .SDcols = cols]

## join data and fix a few FT names
cols <- grep("Stand.age.1", names(patchVegTypeCN), invert = TRUE, value = TRUE)   ## we don't have to bring the Stand.age.1 as it's the same as calculated above
ageDataCN <- patchVegTypeCN[, ..cols][ageDataCN, on = .(Patch.ID)]
ageDataCN$Cover.dendro <- sub("Mixedwood", "mixedwood", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("Broadleaf", "broadleaf", ageDataCN$Cover.dendro)
ageDataCN$Cover.dendro <- sub("-", "", ageDataCN$Cover.dendro)
rm(patchVegTypeCN)

## create LandR species code
ageDataCN[grepl("POBA|POTR", Species), speciesCode := "Popu_sp"]  ## collapse populus
ageDataCN[grepl("ABLA", Species), speciesCode := "Abie_sp"]
ageDataCN[grepl("PIEN", Species), speciesCode := "Pice_eng"]
ageDataCN[grepl("PICO", Species), speciesCode := "Pinu_sp"]
ageDataCN[grepl("PSME", Species), speciesCode := "Pseu_men"]

## calculate no. cohorts:
ageDataCN[, noCohorts := length(unique(Est.bin)) , by = .(Patch.ID)]

## make fire presence/absence column
## all patches with mean.FI == NA had one inferred scar (there are no plots without fires)
ageDataCN[, firePresAbs := 1]

## calculate age at minimum DBH (4cm) -- used to filter simulated data

modelData <- rbind(ageDataCN[, .(speciesCode, DBH.cm, Reconstructed.age)],
                   expand.grid(Reconstructed.age = 1,
                               DBH.cm = 1,
                               speciesCode = na.omit(unique(ageDataCN$speciesCode))))   ## force the intercept at 0 by adding 0 data
modAgeDBH <- lm(Reconstructed.age ~ DBH.cm*speciesCode,
                data = modelData) ## na's are remnants

modAgeDBH2 <- lm(log(Reconstructed.age) ~ DBH.cm*speciesCode,
                 data = modelData) ## na's are remnants

modAgeDBH3 <- lm(log(Reconstructed.age) ~ log(DBH.cm)*speciesCode,
                 data = modelData) ## na's are remnants

if (FALSE) {
  performance::check_model(modAgeDBH)
  performance::check_model(modAgeDBH2)
  performance::check_model(modAgeDBH3)


  AICcmodavg::AICc(modAgeDBH)
  AICcmodavg::AICc(modAgeDBH2)
  AICcmodavg::AICc(modAgeDBH3)  ## best

  plotData <- na.omit(ageDataCN[, .(speciesCode, DBH.cm, Reconstructed.age)])
  plotData <- rbind(plotData, expand.grid(Reconstructed.age = NA,
                                          DBH.cm = c(1, 4),
                                          speciesCode = unique(plotData$speciesCode)))
  plotData[, `:=`(pred = predict(modAgeDBH, plotData),
                  predLog = predict(modAgeDBH2, plotData),
                  predLogLog = predict(modAgeDBH3, plotData))]
  ggplot(plotData,
         aes(x = DBH.cm, y = Reconstructed.age)) +
    geom_point(alpha = 0.5) +
    geom_point(data = data.table(DBH.cm = 0, Reconstructed.age = 0),
               colour = "grey") +
    # stat_smooth(method = "lm", colour = "blue") +
    # geom_line(aes(y = pred), col = "red") +
    # geom_line(aes(y = exp(predLog)), col = "green") +
    geom_line(aes(y = exp(predLogLog)), col = "blue") +
    theme_pubr() +
    facet_wrap(~ speciesCode) +
    labs(x = "DBH (cm)", y = "Reconstructed age")
  ggsave(file.path(figOutputPath, "ageMinDBHmodel.png"), width = 10, height = 7, dpi = 300)

  rm(plotData)
}

## predict age at minimum DBH (4 cm)
newData <- data.table(DBH.cm = 4, speciesCode = unique(ageDataCN[!is.na(speciesCode)]$speciesCode))
newData[, ageAtMinDBH := round(exp(predict(modAgeDBH3, newdata = newData)), 0)]
newData <- newData[complete.cases(newData)]

## check -- are predictions larger than the ages of the smallest trees? no, but only IF 0,0 data is added and log-log model is used
ageDataCN[, Reconstructed.age[which(DBH.cm == min(DBH.cm, na.rm = TRUE))], by = speciesCode]

ageDataCN <- newData[, .(speciesCode, ageAtMinDBH)][ageDataCN, on = .(speciesCode)]

rm(newData)
