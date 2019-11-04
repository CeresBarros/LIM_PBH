## --------------------------------
## MODULE HEX STICKERS
## --------------------------------

library("hexSticker")
library("showtext")
library("sysfonts")
library("reproducible")

if (!file.exists("E:/GitHub/LandscapesInMotion/R/SpaDES/moduleStickers/moduleTable.csv"))
  moduleTable <- prepInputs(targetFile = "moduleTable.csv",
                            url = "https://github.com/tati-micheletti/host/raw/master/stickers/moduleTable.csv",
                            destinationPath = "E:/GitHub/LandscapesInMotion/R/SpaDES/moduleStickers",
                            fun = "read.csv") else {
                              moduleTable <- as.data.table(read.csv("E:/GitHub/LandscapesInMotion/R/SpaDES/moduleStickers/moduleTable.csv",
                                                                    stringsAsFactors = FALSE))
                            }


moduleList <- unlist(modules(simList_noPM))
names(moduleList) <- NULL
moduleList[6] <- "fireSpread"

lapply(moduleList, moduleSticker, moduleTable = moduleTable,
       directory = "E:/GitHub/LandscapesInMotion/R/SpaDES/moduleStickers/",
       imageFolder = "E:/GitHub/LandscapesInMotion/R/SpaDES/moduleStickers/")
