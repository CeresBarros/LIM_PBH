## -----------------------------------
## LOAD/MAKE NECESSARY OBJECTS
## -----------------------------------

## this script loads/treats/makes the necessary objects for the simulation
## it should be sourced before running any modules

## STUDY AREA(S) ---------------------------------------

## Foothills and a smaller region for testing
## prepInputs doens't work with kmz, so download and unzipping need to be done externally.
foothills <- Cache(prepKMZ2shapefile,
                   url = "https://drive.google.com/open?id=1OCqRRIjRNFi6LmxY6m8QH4gMBOLTNeDs",
                   archive = "Foothills_study_area.zip",
                   destinationPath = "data/maps",
                   cacheRepo = "data/cache",
                   userTags = "foothills",
                   omitArgs = c("userTags"))
foothills <- sp::spTransform(foothills,
                             sp::CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
foothillsSMALL <- raster::buffer(foothills, width = -30000)
foothillsMED <- raster::buffer(foothills, width = -15000)

## RASTER TO MATCH LARGE AND KNN 2010 THINGS -----------
## get rawBiomassMap first
rawBiomassMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                           "canada-forests-attributes_attributs-forests-canada/",
                           "2011-attributes_attributs-2011/",
                           "NFI_MODIS250m_2011_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")
standAgeMapURL <- paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                         "canada-forests-attributes_attributs-forests-canada/",
                         "2011-attributes_attributs-2011/",
                         "NFI_MODIS250m_2011_kNN_Structure_Stand_Age_v1.tif")
fireURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"

cacheTags <- c("1_simObjects")

rawBiomassMapFilename <- basename(rawBiomassMapURL)
rawBiomassMap <- Cache(prepInputs,
                       targetFile = rawBiomassMapFilename,
                       url = rawBiomassMapURL,
                       destinationPath = simPaths$inputPath,
                       studyArea = foothills,
                       rasterToMatch = NULL,
                       maskWithRTM = FALSE,
                       useSAcrs = FALSE,     ## never use SA CRS
                       method = "bilinear",
                       datatype = "INT2U",
                       filename2 = NULL,
                       cacheRepo = simPaths$cachePath,
                       userTags = c(cacheTags, "rawBiomassMap"),
                       omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))

## now RTMLarge
rasterToMatchLarge <- rawBiomassMap
RTMvals <- raster::getValues(rasterToMatchLarge)
rasterToMatchLarge[!is.na(RTMvals)] <- 1
rasterToMatchLarge <- Cache(writeOutputs, rasterToMatchLarge,
                            filename2 = file.path(simPaths$cachePath, "rasters", "rasterToMatchLarge.tif"),
                            datatype = "INT2U", overwrite = TRUE,
                            userTags = c(cacheTags, "rasterToMatchLarge"),
                            cacheRepo = simPaths$cachePath,
                            omitArgs = c("userTags"))

## and standAgeMap
standAgeMap <- Cache(LandR::prepInputsStandAgeMap,
                     destinationPath = simPaths$inputPath,
                     ageURL = standAgeMapURL,
                     studyArea = raster::aggregate(foothills),
                     rasterToMatch = rasterToMatchLarge,
                     filename2 = .suffix("standAgeMap.tif", paste0("_", reproducible::studyAreaName(foothills))),
                     overwrite = TRUE,
                     fireURL = fireURL,
                     fireField = "YEAR",
                     startTime = 2011,
                     cacheRepo = simPaths$cachePath,
                     userTags = c("prepInputsStandAge_rtm", cacheTags),
                     omitArgs = c("destinationPath", "targetFile", "overwrite",
                                  "alsoExtract", "userTags"))

## NON FOREST FUELS TABLE LCC2010 ----
## see data/LCC2010_LCC2005_correspondence.xlsx
nonForestFuelsTable <- data.table(LC = c(8, 10, 12, 11, 13, 16),
                                  FuelTypeFBP = c("O1b", "O1b", "O1b", "O1b", "O1b", "NF"),
                                  FuelType = c(15, 15, 15, 15, 15, 19),
                                  fixedCuring = c(FALSE, FALSE, FALSE, FALSE, TRUE, NA),
                                  curingMean = c(60, 60, 60, 60, 30, NA),
                                  curingMin = c(50, 50, 50, 50, 30, NA),
                                  curingMax = c(80, 80, 80, 80, 30, NA))
## ECOREGION LAYER --------------------
ecoregionLayer <- Cache(prepInputs,
                        targetFile = "Natural_Regions_Subregions_of_Alberta.shp",
                        archive = asPath("natural_regions_subregions_of_alberta.zip"),
                        url = "https://www.albertaparks.ca/media/429607/natural_regions_subregions_of_alberta.zip",
                        alsoExtract = "similar",
                        studyArea = foothills,
                        fun = "raster::shapefile",
                        destinationPath = simPaths$inputPath,
                        overwrite = TRUE,
                        useSAcrs = TRUE,
                        cacheRepo = "data/cache",
                        userTags = c("prepInputsNatSubRegionsAB_SA"))
ecoregionLayer$ecozoneCode <- as.numeric(factor(paste(ecoregionLayer$NRNAME, ecoregionLayer$NSRNAME)))

## SPECIES LISTS ------------------------
## Set up sppEquiV
data("sppEquivalencies_CA", package = "LandR")
sppEquivalencies_CA[grep("Pin", LandR), `:=`(EN_generic_short = "Pine",
                                             EN_generic_full = "Pine",
                                             Leading = "Pine leading")]

## Make LIM spp equivalencies column and add correspondences to other columns.
sppEquivalencies_CA[, LIM := c(Abie_bal = "Abie_sp", Abie_las = "Abie_sp", Abie_sp = "Abie_sp",
                               Lari_lar = "Lari_sp", Lari_lya = "Lari_sp",
                               Pice_mar = "Pice_mar", Pice_gla = "Pice_gla", Pice_eng = "Pice_eng",
                               Pinu_alb = "Pinu_sp", Pinu_ban = "Pinu_sp", Pinu_con = "Pinu_sp", Pinu_fle = "Pinu_sp", ## flammable pines
                               Pinu_pon = "Pinu_pon", ## not so flammable pine
                               Popu_tre = "Popu_sp", Betu_pap = "Popu_sp", Popu_bal = "Popu_sp",
                               Pseu_men = "Pseu_men")[LandR]]

sppEquivalencies_CA[, EN_generic_short := c(Abie_sp = "Fir",
                                            Lari_sp = "Tamarack",
                                            Pice_mar = "Bl spruce", Pice_gla = "Wh spruce", Pice_eng = "Eng spruce",
                                            Pinu_sp = "Pine",
                                            Pinu_pon = "Pd pine",
                                            Popu_sp = "Decid",
                                            Pseu_men = "Doug-fir")[LIM]]
sppEquivalencies_CA[, EN_generic_full := c(Abie_sp = "Fir",
                                           Lari_sp = "Tamarack",
                                           Pice_mar = "Black spruce", Pice_gla = "White spruce", Pice_eng = "Engelmann spruce",
                                           Pinu_sp = "Pine",
                                           Pinu_pon = "Ponderosa pine",
                                           Popu_sp = "Deciduous",
                                           Pseu_men = "Douglas-fir")[LIM]]

sppEquivalencies_CA[, FI_layers := c(Abie_sp = "Fir",
                                     Lari_sp = "",
                                     Pice_mar = "Black.Spruce", Pice_gla = "White.Spruce", Pice_eng = "",
                                     Pinu_sp = "Pine",
                                     Pinu_pon = "Pine",
                                     Popu_sp = "Deciduous",
                                     Pseu_men = "")[LIM]]

sppEquivalencies_CA[, Leading := c(Abie_sp = "Fir leading",
                                   Lari_sp = "Tamarack Leading",
                                   Pice_mar = "Black.Spruce", Pice_gla = "White.Spruce", Pice_eng = "",
                                   Pinu_sp = "Pine leading",
                                   Pinu_pon = "Pd pine leadinng",
                                   Popu_sp = "Deciduous leading",
                                   Pseu_men = "Doug-fir leading")[LIM]]

## Make Cameron's spp equivalencies column for species used in LIM.
sppEquivalencies_CA[, Cameron := c(Abie_las = "ABLA", Abie_sp = "ABLA",
                                   Pice_gla = "PIGL", Pice_eng = "PIEN",
                                   Pinu_con = "PICO", Pinu_fle = "PIFL", ## flammable pines
                                   Pinu_pon = "PIPO", ## not so flammable pine
                                   Popu_tre = "POTR", Betu_pap = "BEPA", Popu_bal = "POBA",
                                   Pseu_men = "PSME")[LandR]]

## define spp column to use for model
sppEquivCol <- "LIM"
sppEquivalencies_CA <- na.omit(sppEquivalencies_CA, cols = sppEquivCol)

## PSP DATA ------------------------------------------------------------------
## Set up PSP data for Biomass_speciesParameters
PSPmeasure_sppParams <- Cache(prepInputs,
                              targetFile = "randomizedPSPmeasure.rds",
                              archive = "randomized_Biomass_speciesParameters_Inputs.zip",
                              url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                              fun = "readRDS",
                              destinationPath = simPaths$inputPath,
                              cacheRepo = "data/cache",
                              userTags = "PSPmeasure_sppParams",
                              omitArgs = c("userTags"))

PSPplot_sppParams <- Cache(prepInputs,
                           targetFile = "randomizedPSPplot.rds",
                           archive = "randomized_LandR_speciesParameters_Inputs.zip",
                           url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                           destinationPath = simPaths$inputPath,
                           fun = "readRDS",
                           cacheRepo = "data/cache",
                           userTags = "PSPplot",
                           omitArgs = c("userTags"))

PSPgis_sppParams <- Cache(prepInputs,
                          targetFile = "randomizedPSPgis.rds",
                          archive = "randomized_LandR_speciesParameters_Inputs.zip",
                          url = "https://drive.google.com/file/d/1LmOaEtCZ6EBeIlAm6ttfLqBqQnQu4Ca7/view?usp=sharing",
                          destinationPath = simPaths$inputPath,
                          fun = "readRDS",
                          cacheRepo = "data/cache",
                          userTags = "PSPgis_sppParams",
                          omitArgs = c("userTags"))
