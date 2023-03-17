## GET CAMERON'S AGE DATA AND STAND VEG TYPES -----------------------
source("data/CameronsAgeData/Stand age function_for Ceres_Adapted.R")
ageDataCN <- fread("data/CameronsAgeData/treelist_outputs_for Ceres v2.csv")   ## v2 has stand density and basal area
patchVegTypeCN <- fread("data/CameronsAgeData/patch outputs_for Ceres_Oct2021.csv")

## remove records after 1940 (fire exclusion)
ageDataCN <- ageDataCN[Establishment.date <= 1940]

## remove funky record
ageDataCN <- ageDataCN[Reconstructed.age != 2018]

## calculate stand-ages for FT-level comparisons
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
modAgeDBH <- lm(log(Reconstructed.age) ~ DBH.cm*speciesCode, data = ageDataCN[!is.na(speciesCode)]) ## na's are remnants

if (FALSE) {
  ggplot(ageDataCN[!is.na(speciesCode)],
         aes(x = DBH.cm, y = log(Reconstructed.age), colour = speciesCode)) +
    geom_point(alpha = 0.5) +
    stat_smooth(method = "lm", colour = "blue") +
    theme_pubr() +
    facet_wrap(~ speciesCode)
}

## predict age at minimum DBH (4 cm)
newData <- data.table(DBH.cm = 4, speciesCode = unique(ageDataCN$speciesCode))
newData[, ageAtMinDBH := round(predict(modAgeDBH, newdata = newData), 0)]
ageDataCN <- newData[, .(speciesCode, ageAtMinDBH)][ageDataCN, on = .(speciesCode)]

rm(newData)
