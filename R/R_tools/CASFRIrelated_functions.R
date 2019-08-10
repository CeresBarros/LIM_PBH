## ------------------------------------------------------
## CASFRI-RELATED FUNCTIONS
##
## Ceres: Sep 2018
## ------------------------------------------------------

## NOTE: CASFRI functions follow SK_conversion06MISTIK.pm for SK and AB_conversion31.pm for AbB

## CONVERT VALUES TO CASFIR USING A MAPPING TABLE -------
## dt is a data.table containing the columns to be converted
## correspTab is a correspondence/conversion table
## dtVar is the column contaning the values to be replaced
## correspVar is the column in the correspTab that matches dtVar
## newVar is the the column in correspTab containing the values to use as replacement
## newName is the name to give to the "converted" column, if NULL dtVar will be used
## keepOld determines wheter the original colum is kept or not. Defaults to FALSE, which removes the column

invent2CASFRI <- function(dt, correspTab, dtVar, correspVar, newVar,
                          newName = NULL, keepOld = FALSE) {
  if (is.null(newName)) newName <- dtVar

  ## convert to caps if character
  if (is.character(dt[, ..dtVar])) {
    set(dt, NULL, j = dtVar,
        value = sapply(dt[, ..dtVar], function(x) factor(toupper(x))))
  }
  ## join the new variable by correspondence between
  ## dtVar and correspVar
  onCorresp <- paste0(correspVar, "==", dtVar)
  cols <-  c(newVar, correspVar)
  dt <- correspTab[, ..cols]  %>%
    .[dt, on = onCorresp, nomatch = NA]

  ## rename variables
  setnames(dt, old = newVar, new = newName)
  if (keepOld)
    setnames(dt, old = correspVar, new = dtVar)
  else
    set(dt, NULL, j = correspVar, value = NULL)

  return(dt)
}

## GET SPP LATIN NAMES ---------
## Equivalent to function "Latine" the CASFRI .pm files (AB & SK)
## sp: is the species code, on species at ta time
## spTable: is the table with species name correspondences

spLatinName <- function(sp, spTable,
                        ERRCODE = "XXXX ERRC", MISSCODE = "XXXX MISS") {
  if (!any(grepl("data.table", class(spTable))))
    spTable <- data.table(spTable)
  if (!any(grepl("CASFRI", names(spTable))))
    stop("Can't find a 'CASFRI' column")
  if (ncol(spTable) != 2)
    stop("'spTable' should have 2 columns of spp codes; one has to be named 'CASFRI'")
  if (length(sp) > 1)
    stop("Provide only one species at a time (per layer, per stand)")

  ## svae column name for warning, then change it
  provinceCol <- names(spTable)[!grepl("CASFRI", names(spTable))]
  names(spTable)[grepl(provinceCol, names(spTable))] <- "provinceCol"

  ## convert to CAPS
  sp <- toupper(sp)
  spTable[, provinceCol := toupper(provinceCol)]

  ## check if species exists
  if (!is.na(sp)) {
    if (!any(sp %in% spTable[, provinceCol])) {
      spLatin <- ERRCODE
      warning(paste(sp, "was not found in the", provinceCol, "column of 'spTable'"))
    } else {
      spLatin <- spTable[provinceCol == sp, CASFRI]
    }
  } else
    spLatin <- MISSCODE

  return(as.character(spLatin))
}

## DETERMINE SPECIES TYPE (HARDWOOD/SOFTWOOD)
## Function used to determine whether species is softwood or hardwood
## used for further verification in species percentage determination (comment in SK_conversion06MISTIK.pm)
## TypeForest is the function name in the .pm file
## this function is only used in the UTM type SK forest inventories
## sp: is the species code
## province: full name or acronym of province of the vegetation inventory data

TypeForest <- function(sp, province = NULL) {
  ## checks
  if (length(sp) > 1)
    stop("Provide only one species at a time")
  if (is.null(province))
    stop("Please chose a province (full name or acronym)")

  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("Saskatchewan")))
    stop("Only Saskatchewan has implemented methods so far")

  if (province %in% c("SK", "SASKATCHEWAN")) {
    if(is.na(sp)) {
      spType <- "NA"
    } else {
      sp <- toupper(sp)
      if(grepl("GA|TA|BP|WB|WE|MM|BO", sp))  {
        spType <- "H"
      } else {
        if(grepl("WS|BS|JP|BF|TL|LP", sp)) {
          spType <- "S"
        } else {
          spType <- "NA"
          warning(paste(sp, "is not a recognised Saskatchewan species code"))
        }
      }
    }
  }

  return(spType)

}

## CALCULATE SPECIES PERCENTAGE --------
## spPer: is the original species cover
## province: full name or acronym of province of the vegetation inventory data

spPercent <- function(spPer, province = NULL) {
  ## checks
  if (length(spPer) > 1)
    stop("Provide only one species at a time (per layer, per stand)")
  if (is.null(province))
    stop("Please chose a province (full name or acronym)")

  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("Alberta", "Saskatchewan", "AB", "SK")))
    stop("Only Alberta and Saskatchewan have implemented methods so far")

  if (province %in% c("SK", "SASKATCHEWAN", "AB", "ALBERTA")) {
    if (is.na(spPer))
      spPer <- 0
    spPer * 10
  }
}

##  this funciton adjusts the percentage when <100
## it needs to be applied across all calculated percentages for each stand PER LAYER
## sppPer: vector of percentage cover of all species in a layer
spPercentAdjust <- function(sppPer, province = NULL) {
  ## checks
  sppPer <- unlist(sppPer)
  if (length(sppPer) > 10)
    stop("CASFRI only accepts a maximum of 10 species per layer")
  if (length(sppPer) < 1)
    stop("No species percentage covers supplied")

  if (is.null(province))
    stop("Please chose a province (full name or acronym)")

  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("Alberta", "Saskatchewan", "AB", "SK")))
    stop("Only Alberta and Saskatchewan have implemented methods so far")

  if (province %in% toupper(c("SK", "Saskatchewan"))) {
    if(sum(sppPer[1:5] == 90)) sppPer[6] <- 10

    return(as.list(sppPer))
  }

  if (province %in% toupper(c("AB", "Alberta"))) {
    totalPer <- sum(sppPer)
    if (totalPer == 80 &
        sppPer[1] == 50 & sppPer[2] == 20 & sppPer[3] == 10) {
      sppPer[1] == 60
      sppPer[2] == 30
    } else {
      if (totalPer == 90 &
          sppPer[1] == 90) {
        sppPer[1] == 100
      } else {
        if (totalPer == 90 &
            sppPer[1] == 80 & sppPer[2] == 10) {
          sppPer[1] == 90
        } else {
          if (totalPer == 110 &
              sppPer[1] == 60 & sppPer[2] == 30 & sppPer[3] == 20) {
            sppPer[1] == 50
          } else {
            if (totalPer == 40 &
                sppPer[1] == 0 & sppPer[2] == 20 & sppPer[3] == 20) {
              sppPer[1] == 60
            }
          }
        }
      }
    }

    return(as.list(sppPer))
  }
}


## CALCULATE YEAR OF ORIGIN UPPER/LOWER --------
## this function should be applied to each layer
## year: is the original year

originUpper <- function(year, province = NULL,
                        MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks
  if (length(year) > 1)
    stop("Provide only one year at a time")
  if (is.null(province))
    stop("Please chose a province (full name or acronym)")
  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("Alberta", "Saskatchewan", "AB", "SK")))
    stop("Only Alberta and Saskatchewan have implemented methods so far")

  year <- if (province %in% c("SK", "SASKATCHEWAN")) {
    if (!is.na(year) & year > 0) {
      if(year %% 10 > 0) year else year + 5
    } else if (is.na(year)) MISSCODE else ERRCODE
  } else if (province %in% c("AB", "ALBERTA")) {
    if (!is.na(year) & year > 0) {
      if (nchar(as.character(year)) == 4)
        year
      else if (nchar(as.character(year)) == 2)
        as.numeric(paste0(1, year, 9))
    } else if (is.na(year) | year == 0) MISSCODE else ERRCODE
  }
  return(year)
}

originLower <- function(year,  province = NULL,
                        MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks
  if (length(year) > 1)
    stop("Provide only one year at a time")
  if (is.null(province))
    stop("Please chose a province (full name or acronym)")
  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("Alberta", "Saskatchewan", "AB", "SK")))
    stop("Only Alberta and Saskatchewan have implemented methods so far")

  year <- if (province %in% c("SK", "SASKATCHEWAN")) {
    if (!is.na(year) & year > 0) {
      if(year %% 10 > 0) year else year -4
    } else if (is.na(year)) MISSCODE else ERRCODE
  } else if (province %in% c("AB", "ALBERTA")) {
    if (!is.na(year) & year > 0) {
      if (nchar(as.character(year)) == 4)
        year
      else if (nchar(as.character(year)) == 2)
        as.numeric(paste0(1, year, 0))
    } else if (is.na(year) | year == 0) MISSCODE else ERRCODE
  }
  return(year)
}

## ASSESS NATURALLY NON-VEGETATED CLASS SASKATCHEWAN --------
## this function must be applied per stand not per layer, and only to the first 3 layers
## it follows SK_conversion06.pm and the CASFRI manual
## nonForest: is NVSL in SFVI
## aquatic: is AQUATIC CLASS in SFVI (in Dave's data use TYPE)
##        of the three canopy layers (not CASFRI CC)

nonVegNatSK <- function(nonForest, aquatic,
                        MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(unique(nonForest)) > 1 &
      length(unique(aquatic)) > 1)
    stop("there should only be one nonForest/modifier/aquatic value per stand")

  if (length(nonForest) > 3 |
      length(aquatic) > 3)
    stop("suply values from canopy layers only (up to 3 allowed)")

  reps <- length(nonForest)
  nonForest <- toupper(unique(nonForest))
  aquatic <- toupper(unique(aquatic))

  natNonVeg <- if (is.na(aquatic)) nonForest else aquatic

  if(!is.na(natNonVeg)) {
    ## in SK_conversion06MISTIK.pm, SFVI codes are c("L", "R", "FL") probably due to an error in that inventory
    ## other SK scripts have LA, RI, FL, RK staying the same
    if (!natNonVeg %in% c("LA", "RI", "FL", "RK")) {
      if (natNonVeg %in% c("L", "R", "SF", "CB", "MS", "SB", "SA",
                           "UK", "GR", "WA", "ST", "SL", "FP")) {
        corresp <- matrix(c("LA", "RI", "FL", "SL", "EX", "RF", "BE",
                            "OT", "WS", "LA", "RI", "FL", "BP"), nrow = 13, ncol = 1,
                          dimnames = list(c("L", "R", "SF", "CB", "MS", "SB", "SA",
                                            "UK", "GR", "WA", "ST", "SL", "FP")))
        natNonVeg <- corresp[natNonVeg,]
        names(natNonVeg) <- NULL
      } else
        natNonVeg <- ERRCODE
    }
  } else
    natNonVeg <- MISSCODE

  as.character(rep(natNonVeg, reps))
}

## ASSESS ANTHROPOGENIC NON-VEGETATED CLASS SASKATCHEWAN--------
## this function must be applied per stand not per layer, across all layers (i.e. canopy shrub and herb)
## landUse: is LUC in SFVI
## nonForest: is NVSL in SFVI

nonVegAnthSK <- function(landUse, nonForest,
                         MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(unique(landUse)) > 1 &
      length(unique(nonForest)) > 1)
    stop("there should only be one landUse/nonForest unique value per stand")

  if (length(landUse) > 5 |
      length(nonForest) > 5)
    stop("there should only be 5 different veg. layers")

  reps <- length(nonForest)
  landUse[landUse %in% 0] <- NA
  landUse <- toupper(as.character((unique(landUse))))
  nonForest <- toupper(as.character((unique(landUse))))

  landUseCodes <- c("GSOF", "MUOU", "MG", "LMBY",
                    "SGDU", "PEATC", "BUPO", "RWGU",
                    "VEGU", "BUPG", "CMTY", "DMGU",
                    "FTOW")
  if (!is.na(landUse)) {
    if (landUse %in% landUseCodes) {
      anthNonVeg <- landUse
      corresp <- matrix(c(rep("IN", 7), "FA", "CL",
                          "SE", "SE", "OT", "OT"), nrow = 13, ncol = 1,
                        dimnames = list(landUseCodes))
      anthNonVeg <- corresp[anthNonVeg,]
      names(anthNonVeg) <- NULL
    } else
      anthNonVeg <- ERRCODE
  } else if (!is.na(nonForest)) {
    anthNonVeg <- "FA"
  }

  if (!exists("anthNonVeg"))
    anthNonVeg <- MISSCODE

  as.character(rep(anthNonVeg, reps))

}

## ASSESS VEGETATED NON-FORESTED CLASS SASKATCHEWAN --------
## this function must be applied per stand not per layer, across all layers (i.e. canopy shrub and herb)
## sp1: is the dominant species, only species from shrub and herb layers are considered
## cover: is a vector of original crown closure across all layers (not CASFRI CC)
## moist: is a vector of original soil moisture regime cross all layers
## layerID: is a vector of layer type indices. In Saskatchewan 1:3 are crown layers 1, 2 and 3,
##      4 is the shrub layer and 5 is the herbaceous layer (LAYER_RANK is casfri).

nonForestVegSK <- function(sp1, cover, moist, layerID,
                           MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(sp1) > 5 |
      length(cover) > 5 |
      length(moist) > 5)
    stop("there should only be 5 different layers")

  if (length(unique(moist)) > 1)
    stop("there should only be one moist unique value per stand")

  moist <- toupper(as.character(unique(moist)))

  ## get shurb and herb leading species
  sp1Shrub <- toupper(sp1[layerID == 4])
  sp1Herb <- toupper(sp1[layerID == 5])

  ## make separate cover objects for each layer
  cover[is.na(cover)] <- 0
  cover1 <- cover[layerID == 1]
  cover2 <- cover[layerID == 2]
  cover3 <- cover[layerID == 3]
  coverShrub <- cover[layerID == 4]
  coverHerb <- cover[layerID == 5]

  ## calculate the sum of canopy cover
  sumCC <- sum(cover1, cover2, cover3)

  if(sumCC > 0 & sumCC < 10) {
    if(sum(coverShrub, coverHerb) > 0) {
      if (coverShrub >= coverHerb) {
        nonForVeg <- sp1Shrub
        if(!is.na(nonForVeg)) {
          corresp <- matrix(c(rep("ST", 8), rep("SL", 16)), nrow = 24, ncol = 1,
                            dimnames = list(c("TS", "AL", "BH", "MA",
                                              "SA", "PC", "CR", "WI",
                                              "LS", "RO", "BI", "BU",
                                              "DW", "RA", "CU", "SN",
                                              "BB", "CI", "BL", "LA",
                                              "LE", "BE", "LC", 'LB')))
          nonForVeg <- corresp[nonForVeg,]
          names(nonForVeg) <-  NULL
        }  else nonForVeg <- MISSCODE
      } else if (coverShrub < coverHerb) {
        nonForVeg <- sp1Herb
        if(!is.na(nonForVeg)) {
          corresp <- matrix(c(rep("HF", 7), c("HE", "HG", "HG", "BR", "BR", "BR")),
                            nrow = 13, ncol = 1,
                            dimnames = list(c("FE", "SA", "BU", "ST",
                                              "NC", "MM", "CO",
                                              "HE", "GR", "SE", "FE",
                                              "SM", "LI")))
          nonForVeg <- corresp[nonForVeg,]
          names(nonForVeg) <-  NULL
        }  else nonForVeg <- MISSCODE
      }
    } else {
      nonForVeg <- moist
      if (!is.na(nonForVeg)) {
        corresp <- matrix(rep("OM", 4), nrow = 4, ncol = 1,
                          dimnames = list(c("VM", "MW", "W", "VW")))
        nonForVeg <- corresp[nonForVeg,]
        names(nonForVeg) <-  NULL
      }  else nonForVeg <- MISSCODE
    }
  } else nonForVeg <- ERRCODE

  as.character(rep(nonForVeg, length(cover)))

}

## WETLAND CODES SASKATCHEWAN --------
## this function must be applied per stand not per layer, across all layers (i.e. canopy shrub and herb)
## moist: is a vector of original soil moisture regime cross all layers (not CASFRI SMR)
## cover: is a vector of  original crown closure (although only first layer is used)  (not CASFRI CC)
## height: is a vector of original height (although only first layer is used) (not CASFRI CC)
## nonForest: is a vector of NVSL in SFVI (should only have one unique value)
## sp1: is a vector of original dominant species (although only first layer is used) (original data, not CAFRI SPEC1)
## sp2: is a vector of the second most-dominant species (although only first layer is used)  (original data, not CAFRI SPEC2)
## sp1Per: is a vector of the % cover of sp1, already calculated using the CASFRI approach (although only first layer is used) (=SPEC1_PER in layer 1)
## layerID: is a vector of layer type indices. In Saskatchewan 1:3 are crown layers 1, 2 and 3,
##      4 is the shrub layer and 5 is the herbaceous layer (LAYER_RANK is casfri).

wetlandCodesSK <- function(moist, cover, height, nonForest, sp1,
                           sp2, sp1Per,  layerID, province = NULL,
                           MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(sp1) > 5 |
      length(cover) > 5 |
      length(moist) > 5 |
      length(sp2) > 5 |
      length(nonForest) > 5 |
      length(sp1Per) > 5 |
      length(layerID) > 5)
    stop("Saskatchewan: there should only be 5 different veg. layers")

  moist <- toupper(as.character(moist[layerID == 1]))
  nonForest <- toupper(as.character(nonForest[layerID == 1]))
  sp1 <- toupper(as.character(sp1[layerID == 1]))
  sp2 <- toupper(as.character(sp2[layerID == 1]))
  cover <- as.numeric(cover[layerID == 1])
  height <- as.numeric(height[layerID == 1])
  sp1Per <- as.numeric(sp1Per[layerID == 1])

  wetCodes <- list(WETLAND_CLASS = NA, WETLAND_VEG_MOD = NA,
                   WETLAND_LAND_MOD = NA, WETLAND_LOCAL_MOD = NA)
  if (is.na(sp1Per))
    sp1Per <- 0

  if (moist %in% c("MW", "W", "VW")) {
    if(nonForest %in% c("HE", "GR", "MO", "AV", "TS", "LS")) {
      corresp <- matrix(c("M,O,N,G", "M,O,N,G", "F,O,N,G", "O,O,N,N", "S,O,N,S", "S,O,N,S"),
                        nrow = 6, ncol = 1,
                        dimnames = list(c("HE", "GR", "MO", "AV", "TS", "LS")))
      wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[nonForest,], split = ",")[[1]])
    } else {
      if (moist %in% "MW") {
        if (sp1 %in% "BS" & sp1Per == 100 &
            cover <= 50 & height < 12) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("B,T,N,N", split = ",")[[1]])
        } else {
          wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
        }
      }

      if (moist %in% "W" & sp1 %in% "BS" & sp1Per == 100) {
        if (cover <= 50 & height < 12) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("B,T,N,N", split = ",")[[1]])
        }
        if (cover < 70 & height >= 12) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
        }
        if (cover >= 70 & height >= 12) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
        }
      }

      if (moist %in% "VW") {
        if (sp1 %in% c("BS", "TL", "WB", "MM", "BP") &
            sp2 %in%  c("BS", "TL", "WB", "MM", "BP") & cover >= 50) {
          if (cover < 70 & height >= 12) {
            wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
          }
          if (cover >= 70) {
            wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
          }
        }
        if (sp1 %in% c("BS", "TL") &
            sp2 %in%  c("BS", "TL") &
            cover <= 50) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("F,T,N,N", split = ",")[[1]])
        } else {
          if (sp1 %in% "TL" & sp1Per == 100) {
            if (cover > 50 & cover < 70 & height >= 12) {
              wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
            } else if (cover >= 70) {
              wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
            } else if (cover <= 50) {
              wetCodes[names(wetCodes)] <- as.list(strsplit("F,T,N,N", split = ",")[[1]])
            }
          }
        }

        if (sp1 %in% c("GA", "WE", "WB", "MM") & sp1Per == 100) {
          if (cover < 70) {

            wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])

          } else {
            wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
          }
        }
      }
    }
  }

  ## if no conditions were met for wetlands attribute missing code
  if (all(sapply(wetCodes, is.na))) {
    wetCodes <- lapply(wetCodes, FUN = function(x) as.character(MISSCODE))
  }

  ## repeat values for all layers
  wetCodes <- lapply(wetCodes, FUN = function(x) rep(x, length(layerID)))
  return(wetCodes)
}

## because Dave's data has no info on NVSL, use stand type instead
wetlandCodesSK2 <- function(moist, cover, height, standType, sp1,
                            sp2, sp1Per,  layerID, province = NULL,
                            MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(sp1) > 5 |
      length(cover) > 5 |
      length(moist) > 5 |
      length(sp2) > 5 |
      length(standType) > 5 |
      length(sp1Per) > 5 |
      length(layerID) > 5)
    stop("Saskatchewan: there should only be 5 different veg. layers")

  if(length(unique(standType)) > 1)
    stop("Saskatchewan: there should only be one stand type value per stand")
  if(length(unique(moist)) > 1)
    stop("Saskatchewan: there should only be one moisture value per stand")

  moist <- toupper(as.character(moist[layerID == 1]))
  standType <- toupper(as.character(standType[layerID == 1]))
  sp1 <- toupper(as.character(sp1[layerID == 1]))
  sp2 <- toupper(as.character(sp2[layerID == 1]))
  cover <- as.numeric(cover[layerID == 1])
  height <- as.numeric(height[layerID == 1])
  sp1Per <- as.numeric(sp1Per[layerID == 1])

  wetCodes <- list(WETLAND_CLASS = NA, WETLAND_VEG_MOD = NA,
                   WETLAND_LAND_MOD = NA, WETLAND_LOCAL_MOD = NA)
  if (is.na(sp1Per))
    sp1Per <- 0

  ## adapted from SFVI for stand type (see wetlandCodesSK and SK_conversion06MISTIK.pm)
  if (moist %in% c("MW", "W", "VW") &
      standType %in% c("WAT", "TMS", "BSH", "OMS", "GRS")) {
    corresp <- matrix(c("O,O,N,N", "W,T,-,-", "S,O,N,S", "M,O,N,G", "M,O,N,G"),
                      nrow = 5, ncol = 1,
                      dimnames = list(c("WAT", "TMS", "BSH", "OMS", "GRS")))
    wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[standType,], split = ",")[[1]])
  }

  if (moist %in% "MW") {
    if (sp1 %in% "BS" & sp1Per == 100 &
        cover <= 50 & height < 12) {
      wetCodes[names(wetCodes)] <- as.list(strsplit("B,T,N,N", split = ",")[[1]])
    } else {
      wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
    }
  }

  if (moist %in% "W" & sp1 %in% "BS" & sp1Per == 100) {
    if (cover <= 50 & height < 12) {
      wetCodes[names(wetCodes)] <- as.list(strsplit("B,T,N,N", split = ",")[[1]])
    }
    if (height >= 12) {
      wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
    }
  }


  if (moist %in% "VW") {
    if (sp1 %in% c("BS", "TL", "WB", "MM", "BP") &
        sp2 %in%  c("BS", "TL", "WB", "MM", "BP") & cover >= 50) {
      if (cover < 70 & height >= 12) {
        wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
      }
      if (cover >= 70) {
        wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
      }
    }
    if (sp1 %in% c("BS", "TL") &
        sp2 %in%  c("BS", "TL") &
        cover <= 50) {
      wetCodes[names(wetCodes)] <- as.list(strsplit("F,T,N,N", split = ",")[[1]])
    } else {
      if (sp1 %in% "TL" & sp1Per == 100) {
        if (cover > 50 & cover < 70 & height >= 12) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("S,T,N,N", split = ",")[[1]])
        } else if (cover >= 70) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
        } else if (cover <= 50) {
          wetCodes[names(wetCodes)] <- as.list(strsplit("F,T,N,N", split = ",")[[1]])
        }
      }
    }

    if (sp1 %in% c("GA", "WE", "WB", "MM") & sp1Per == 100) {
      wetCodes[names(wetCodes)] <- if (cover < 70) {
        as.list(strsplit("S,T,N,N", split = ",")[[1]])
      } else {
        wetCodes[names(wetCodes)] <- as.list(strsplit("S,F,N,N", split = ",")[[1]])
      }
    }
  }

  ## if no conditions were met for wetlands attribute missing code
  if (all(sapply(wetCodes, is.na))) {
    wetCodes <- lapply(wetCodes, FUN = function(x) as.character(MISSCODE))
  }

  ## repeat values for all layers
  wetCodes <- lapply(wetCodes, FUN = function(x) rep(x, length(layerID)))
  return(wetCodes)
}

## WETLAND CODES ALBERTA --------
## this function must be applied per stand and accross all layers although only the first layer values are used
## moist: is a vector of CASFRI SMR codes across layers
## cover: is a vector of is the original crown closure of vegetation layers (although only first layer is used) (not CASFRI CC)
## nonForestLand: is a vector of  NFL in AVI across layers
## natNonVeg: is a vector of  NAT_NON in AVI across layers
## sp1: is a vector of dominant species from canopy layer (although only first layer is used)  (original data, not CASFRI SPEC1)
## sp2: is a vector of most-dominant species from canopy layer (although only first layer is used)  (original data, not CASFRI SPEC2)
## sp1Per:  is a vector of  % cover of the sp1 from canopy layer (although only first layer is used)  (original data not CASFRI SPEC1_PER)
## layerID: is a vector of layer type indices. In Alberta 1 is canopy 2 is understory

wetlandCodesAB <- function(moist, cover, nonForestLand, natNonVeg,
                           sp1, sp2, sp1Per, layerID,
                           MISSCODE = -1111L, ERRCODE = -9999L) {
  ## checks:
  if (length(moist) > 2 |
      length(cover) > 2 |
      length(nonForestLand) > 2 |
      length(natNonVeg) > 2 |
      length(sp1) > 2 |
      length(sp2) > 2 |
      length(sp1Per) > 2 |
      length(layerID) > 2)
    stop("there should only be 2 different veg. layers")

  moist <- toupper(as.character(moist[layerID == 1]))
  cover <- toupper(as.character(cover[layerID == 1]))
  nonForestLand <- toupper(as.character(nonForestLand[layerID == 1]))
  natNonVeg <- toupper(as.character(natNonVeg[layerID == 1]))
  sp1 <- toupper(as.character(sp1[layerID == 1]))
  sp2 <- toupper(as.character(sp2[layerID == 1]))
  sp1Per <- as.numeric(sp1Per[layerID == 1])

  wetCodes <- list(WETLAND_CLASS = NA, WETLAND_VEG_MOD = NA,
                   WETLAND_LAND_MOD = NA, WETLAND_LOCAL_MOD = NA)

  if (is.na(sp1Per))
    sp1Per <- 0

    if (moist %in% "W") {
      if (is.na(nonForestLand)) {
        wetCodes <- lapply(wetCodes, FUN = function(x) as.character(MISSCODE))
      } else {
        if (nonForestLand %in% c("S,O,N,S", "S,O,N,S", "M,O,N,G", "M,O,N,G", "F,O,N,G")) {
          corresp <- matrix(c("S,O,N,S", "S,O,N,S", "M,O,N,G", "M,O,N,G", "F,O,N,G"),
                            nrow = 5, ncol = 1,
                            dimnames = list(c("SO", "SC", "HG", "HF", "BR")))
          wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[nonForestLand,], split = ",")[[1]])
        } else {
          if (natNonVeg %in% "NMB")
            wetCodes[names(wetCodes)] <- as.list(strsplit("S,O,N,S", split = ",")[[1]])
        }
      }
    } else {
      if ("LT" %in% c(sp1, sp2)) {
        corresp <- matrix(c("F,T,N,N", "F,T,N,N", "S,T,N,N", "S,F,N,N"),
                          nrow = 4, ncol = 1,
                          dimnames = list(c("A", "B", "C", "D")))
        wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[cover,], split = ",")[[1]])
      } else {
        if (sp1 %in% "SB" && sp1Per == 100) {
          corresp <- matrix(c("B,T,N,N", "B,T,N,N", "S,T,N,N", "S,F,N,N"),
                            nrow = 4, ncol = 1,
                            dimnames = list(c("A", "B", "C", "D")))
          wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[cover,], split = ",")[[1]])
        } else {
          if ((sp1 %in% c("SB", "FB") & !sp2 %in% "LT") |
              sp1 %in% c("SW", "BW", "PB")) {
            corresp <- matrix(c("S,T,N,N", "S,T,N,N", "S,T,N,N", "S,F,N,N"),
                              nrow = 4, ncol = 1,
                              dimnames = list(c("A", "B", "C", "D")))
            wetCodes[names(wetCodes)] <- as.list(strsplit(corresp[cover,], split = ",")[[1]])
          }
        }
      }
    }

  ## if no conditions were met for wetlands attribute missing code
  if (all(sapply(wetCodes, is.na))) {
    wetCodes <- lapply(wetCodes, FUN = function(x) as.character(MISSCODE))
  }

  ## repeat values for all layers
  wetCodes <- lapply(wetCodes, FUN = function(x) rep(x, length(layerID)))
  return(wetCodes)
}


## SOIL MOISTURE REGIME ---------
## For some provinces soil moisture regime is not a mere conversion of codes
## sp1: is the dominant species, only species from the canopy layer are considered (use original or casfri names)
## moist: is a vector of original soil moisture regime cross all layers

soilMoistureRegime <- function(moist, sp1, province = NULL,
                               MISSCODE = -1111L, ERRCODE = -9999L) {

  if (is.null(province))
    stop("Please chose a province (full name or acronym)")

  province <- toupper(province)

  if (!province %in% toupper(c("Alberta", "British Columbia",
                               "Manitoba", "New Brunswick",
                               "Newfoundland and Labrador", "Northwest Territories",
                               "Nova Scotia", "Nunavut",
                               "Ontario", "Prince Edward Island",
                               "Quebec", "Saskatchewan", "Yukon",
                               "AB", "BC", "MB", "NB",
                               "NL", "NT", "NS", "NU",
                               "ON", "PE", "QC", "SK", "YT")))
    stop("province is not a recognized Canadian province")

  if (!province %in% toupper(c("ALBERTA", "AB")))
    stop("Only Alberta has implemented methods so far")

  if (province %in% c("AB", "ALBERTA")) {

    ## checks:
    if (length(sp1) > 1 |
        length(moist) > 1)
      stop("Alberta: please provide sp1 and moist values one layer at a time")

    sp1 <- toupper(as.character(sp1))
    moist <- toupper(as.character(unique(moist)))

    if (is.na(moist)) {
      moistCode <- as.character(MISSCODE)
    } else {
      if (moist %in% c("N", "U")) {
        moistCode <- if (sp1 %in% c("LT", "LARI LARI")) {
          "W"
        } else "M"
      } else {
        corresp <- matrix(c("-1", "D", "F", "W", "A"), nrow = 5, ncol = 1,
                          dimnames = list(c("0", "D", "M", "W", "A")))
        moistCode <- corresp[moist,]
        names(moistCode) <- NULL
      }
    }

    if(!exists("moistCode"))
      moistCode <- as.character(ERRCODE)
  }
  moistCode
}


## ADJUST NFL ALBERTA --------
## this function must be applied per stand and accross all layers
## nonForestLand: is a vector of  NFL in AVI across layers
## sp1: is a vector of dominant species from canopy layer (although only first layer is used)  (original data, not CASFRI SPEC1)
## layerID: is a vector of layer type indices. In Alberta 1 is canopy 2 is understory

NFLAdjustAB <- function(nonForestLand, sp1, layerID) {
  ## checks:
  if (length(nonForestLand) > 2 |
      length(sp1) > 2 |
      length(layerID) > 2)
    stop("there should only be 2 different veg. layers")

  if (sp1[layerID == 1] %in% c("SC", "SO") &
      is.na(nonForestLand[layerID == 1]))
    nonForestLand[layerID == 2] <- sp1[layerID == 1]

  if (sp1[layerID == 2] %in% c("SC", "SO") &
      is.na(nonForestLand[layerID == 2]))
    nonForestLand[layerID == 2] <- sp1[layerID == 2]

  return(nonForestLand)
}


## ADJUST SMR ALBERTA --------
## this function must be applied per stand and accross all layers
## moist: is a vector of CASFRI SMR codes across layers
## layerID: is a vector of layer type indices. In Alberta 1 is canopy 2 is understory

SMRAdjustAB <- function(moist, layerID) {
  ## checks:
  if (length(moist) > 2 |
      length(layerID) > 2)
    stop("there should only be 2 different veg. layers")

  if (moist[layerID == 2] %in% "-1")
    moist[layerID == 2] <- moist[layerID == 1]

  return(moist)
}


## CONVERT VEGETATION INVENTORIES FUNCTIONS ----
## inv is a shapefile containing inventory polygons
## tablesDir is the directory containing the conversion tables .xlsx file
## folder is the directory path where the converted shapefiles will be saved
## dim is used for faster caching - omit "inv"

## Alberta
ABToCASFRI <- function(inv, tablesDir, folder, dim) {
  message("Converting Alberta inventory data to CASFRI standard...")
  ## CASFRI VARIABLES
  casfriVars <- c("CROWN_CLOSURE_LOWER", "CROWN_CLOSURE_UPPER", "DIST1",
                  "DIST1_EXTENT_LOWER", "DIST1_EXTENT_UPPER", "DIST1_YEAR",
                  "DIST2", "DIST2_EXTENT_LOWER", "DIST2_EXTENT_UPPER",
                  "DIST2_YEAR", "DIST3", "DIST3_EXTENT_LOWER",
                  "DIST3_EXTENT_UPPER", "DIST3_YEAR", "ECOSITE",
                  "HEIGHT_LOWER", "HEIGHT_UPPER", "LAYER", "LAYER_RANK",
                  "MAPSHEET_ID", "NATURALLY_NON_VEG", "NON_FORESTED_ANTHRO",
                  "NON_FORESTED_VEG", "NUMBER_OF_LAYERS", "ORIGIN_LOWER",
                  "ORIGIN_UPPER", "ORIGINAL_STAND_ID", "PHOTO_YEAR_MAX",
                  "PHOTO_YEAR_MIN", "POLYGON_AREA", "POLYGON_PERIMITER",
                  "PRODUCTIVE_FOREST", "SITE_CLASS", "SITE_INDEX", "SMR",
                  "SPEC1", "SPEC1_PER",
                  "SPEC2", "SPEC2_PER", "SPEC3",
                  "SPEC3_PER", "SPEC4", "SPEC4_PER",
                  "SPEC5", "SPEC5_PER", "SPEC6",
                  "SPEC6_PER", "SPEC7", "SPEC7_PER", "SPEC8",
                  "SPEC8_PER","SPEC9", "SPEC9_PER","SPEC10",
                  "SPEC10_PER", "STAND_STRUCTURE", "STAND_STRUCTURE_PER",
                  "STAND_STRUCTURE_RANGE", "UNPRODUCTIVE_FOREST", "WETLAND_CLASS",
                  "WETLAND_LAND_MOD", "WETLAND_LOCAL_MOD", "WETLAND_VEG_MOD")

  illegalDisturbances <- c("WA", "CS","SA", "TM",
                           "OM", "DS", "AD", "UP",
                           "TL", "PI", "OR", "NC",
                           "ID", "DE", "CY", "CW",
                           "CO", "BD", "AS", "BR",
                           "SU", "SL", "PR", "OC",
                           "MT", "FT", "CL")

  ## CASFRI ERROR/MISSING CODES
  # MISSCODE <- -1111L
  # ERRCODE <- -9999L
  # SPECIES_MISSCODE <- "XXXX MISS"
  # SPECIES_ERRCODE <- "XXXX ERRC"
  # UNDEF <- -8888L

  MISSCODE <- NA
  ERRCODE <- NA
  SPECIES_MISSCODE <- "NA"
  SPECIES_ERRCODE <- "NA"
  UNDEF <- NA

  ## get all necessary conversion tables
  SClassTable <- read.xlsx(tablesDir, sheetName = "SiteClass", header = TRUE) %>%
    data.table(.)
  SStrucTable <- read.xlsx(tablesDir, sheetName = "StandStructTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, AVI)] %>%
    na.omit(.)
  CCUpperTable <- read.xlsx(tablesDir, sheetName = "CrownClosUpperTable", header = TRUE) %>%
    data.table(.)
  CCLowerTable <- read.xlsx(tablesDir, sheetName = "CrownClosLowerTable", header = TRUE) %>%
    data.table(.)
  SppTable <- read.xlsx(tablesDir, sheetName = "SpeciesTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, AVI)] %>%
    na.omit(.)
  NatNVTable <- read.xlsx(tablesDir, sheetName = "NatNonVegTable", header = TRUE) %>%
    data.table(.)
  NonFATable <- read.xlsx(tablesDir, sheetName = "NonForestAnthTable", header = TRUE) %>%
    data.table(.)
  NonFVTable <- read.xlsx(tablesDir, sheetName = "NonForestVegTable", header = TRUE) %>%
    data.table(.)
  distTable <- read.xlsx(tablesDir, sheetName = "DisturbanceTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, AVI)] %>%
    na.omit(.)
  distExtUTable <- read.xlsx(tablesDir, sheetName = "DistExtUpperTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, AVI)] %>%
    na.omit(.)
  distExtLTable <- read.xlsx(tablesDir, sheetName = "DistExtLowerTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, AVI)] %>%
    na.omit(.)

  ## get inventory data
  casfriDT <- st_set_geometry(inv, NULL) %>%
    data.table(.)

  ## LAYER - not necessary but follows Saskatchewan code for consistency
  casfriDT[!is.na(LAYER),
           LAYER := as.character(as.numeric(factor(LAYER, levels = c("1", "2"),
                                                   labels =  c("1", "2")))),
           by = P_ID]
  casfriDT$LAYER <- as.numeric(casfriDT$LAYER) ## can't coerce directly to numeric above.

  ## LAYER RANK
  ## 1 crown layers, 2 for understory
  casfriDT[, LAYER_RANK := LAYER]

  ## SOIL MOISTURE REGIME
  casfriDT[, SMR := soilMoistureRegime(moist = MOIST_REG, sp1 = SP1,
                                       province = "AB", MISSCODE = MISSCODE, ERRCODE = ERRCODE),
           by = c("P_ID", "LAYER")]
  casfriDT[, SMR := SMRAdjustAB(moist = SMR, layerID = LAYER_RANK),
           by = "P_ID"]


  ## SITE CLASS
  casfriDT[, TPR := toupper(TPR)]
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = SClassTable,
                            dtVar = "TPR", correspVar = "AVI",
                            newVar = "CASFRI", newName = "SITE_CLASS", keepOld = TRUE)
  casfriDT[!is.na(TPR) & is.na(SITE_CLASS), SITE_CLASS := as.character(ERRCODE)]
  casfriDT[is.na(SITE_CLASS), SITE_CLASS := as.character(MISSCODE)]


  ## STAND STRUCTURE
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = SStrucTable,
                            dtVar = "STRUC", correspVar = "AVI",
                            newVar = "CASFRI", newName = "STAND_STRUCTURE", keepOld = TRUE)
  casfriDT[!is.na(STRUC) & is.na(STAND_STRUCTURE), STAND_STRUCTURE := as.character(ERRCODE)]
  casfriDT[is.na(STAND_STRUCTURE), STAND_STRUCTURE := as.character(MISSCODE)]

  ## STAND STRUCTURE PERCENT
  casfriDT[STRUC_VAL < 1 | STRUC_VAL > 9 | is.na(STRUC_VAL), STAND_STRUCTURE_PER := 0]
  casfriDT[STRUC_VAL >= 1 & STRUC_VAL <= 9, STAND_STRUCTURE_PER := STRUC_VAL]

  ## CANOPY CLOSURE UPPER & LOWER
  ## UPPER
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = CCUpperTable,
                            dtVar = "DENSITY", correspVar = "AVI",
                            newVar = "CASFRI", newName = "CROWN_CLOSURE_UPPER", keepOld = TRUE)
  casfriDT[!is.na(DENSITY) & is.na(CROWN_CLOSURE_UPPER), CROWN_CLOSURE_UPPER := ERRCODE]
  casfriDT[is.na(CROWN_CLOSURE_UPPER), CROWN_CLOSURE_UPPER := MISSCODE]

  ## LOWER
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = CCLowerTable,
                            dtVar = "DENSITY", correspVar = "AVI",
                            newVar = "CASFRI", newName = "CROWN_CLOSURE_LOWER", keepOld = TRUE)
  casfriDT[!is.na(DENSITY) & is.na(CROWN_CLOSURE_LOWER), CROWN_CLOSURE_LOWER := ERRCODE]
  casfriDT[is.na(CROWN_CLOSURE_LOWER), CROWN_CLOSURE_LOWER := MISSCODE]

  ## STAND HEIGHT UPPER & LOWER
  casfriDT[HEIGHT > 0 & HEIGHT <= 50, `:=` (HEIGHT_UPPER = as.double(HEIGHT),
                                            HEIGHT_LOWER = as.double(HEIGHT))]
  casfriDT[HEIGHT <= 0 | HEIGHT > 50, `:=` (HEIGHT_UPPER = MISSCODE,
                                            HEIGHT_LOWER = MISSCODE)]
  casfriDT[HEIGHT_UPPER > 0, HEIGHT_UPPER := HEIGHT_UPPER + 0.5]
  casfriDT[HEIGHT_LOWER > 0.5, HEIGHT_LOWER := HEIGHT_LOWER - 0.5]
  casfriDT[is.na(HEIGHT), `:=` (HEIGHT_UPPER = MISSCODE,
                                HEIGHT_LOWER = MISSCODE)]
  casfriDT[HEIGHT_UPPER %in% MISSCODE & (is.na(SP1) | SP1 %in% c("SC", "SO")),
           `:=` (HEIGHT_UPPER = UNDEF,
                 HEIGHT_LOWER = UNDEF)]
  casfriDT[HEIGHT_UPPER %in% MISSCODE & !(is.na(SP1) | SP1 %in% c("SC", "SO")) &
             (MOD1 %in% c("BU", "CC", "SN") & MOD1_EXT %in% c(4,5)),
           `:=` (HEIGHT_UPPER = 1,
                 HEIGHT_LOWER = 0)]

  ## SPECIES
  casfriDT[, `:=` (SPEC1 = spLatinName(SP1, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                   SPEC2 = spLatinName(SP2, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                   SPEC3 = spLatinName(SP3, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                   SPEC4 = spLatinName(SP4, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                   SPEC5 = spLatinName(SP5, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                   SPEC6 = as.character(UNDEF),
                   SPEC7 = as.character(UNDEF),
                   SPEC8 = as.character(UNDEF),
                   SPEC9 = as.character(UNDEF),
                   SPEC10 = as.character(UNDEF)),
           by = seq_len(nrow(casfriDT))]

  ## SPECIES PERCENT - CANOPY
  casfriDT[, `:=` (SPEC1_PER = spPercent(SP1_PER, province = "AB"),
                   SPEC2_PER = spPercent(SP2_PER, province = "AB"),
                   SPEC3_PER = spPercent(SP3_PER, province = "AB"),
                   SPEC4_PER = spPercent(SP4_PER, province = "AB"),
                   SPEC5_PER = spPercent(SP5_PER, province = "AB"),
                   SPEC6_PER = UNDEF,
                   SPEC7_PER = UNDEF,
                   SPEC8_PER = UNDEF,
                   SPEC9_PER = UNDEF,
                   SPEC10_PER = UNDEF),
           by = seq_len(nrow(casfriDT))]

  casfriDT[, c("SPEC1_PER", "SPEC2_PER", "SPEC3_PER",
               "SPEC4_PER", "SPEC5_PER") := spPercentAdjust(sppPer = c(SPEC1_PER,
                                                                       SPEC2_PER,
                                                                       SPEC3_PER,
                                                                       SPEC4_PER,
                                                                       SPEC5_PER),
                                                            province = "AB"),
           by = seq_len(nrow(casfriDT))]

  ## ORIGIN UPPER & LOWER
  casfriDT[,  `:=` (ORIGIN_UPPER = as.integer(originUpper(ORIGIN, province= "AB", MISSCODE = MISSCODE, ERRCODE = ERRCODE)),
                    ORIGIN_LOWER = as.integer(originLower(ORIGIN, province= "AB", MISSCODE = MISSCODE, ERRCODE = ERRCODE))),
           by = 1:nrow(casfriDT)]
  ## add error codes for weird years
  casfriDT[(ORIGIN_UPPER > 0 & ORIGIN_UPPER < 1600) | ORIGIN_UPPER > 2014,
           `:=` (ORIGIN_UPPER = ERRCODE, ORIGIN_LOWER = ERRCODE)]


  ## NATURALLY_NON_VEG
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = NatNVTable,
                            dtVar = "NAT_NON", correspVar = "AVI",
                            newVar = "CASFRI", newName = "NATURALLY_NON_VEG", keepOld = TRUE)

  casfriDT[!is.na(NAT_NON) & is.na(NATURALLY_NON_VEG), NATURALLY_NON_VEG := as.character(ERRCODE)]
  casfriDT[is.na(NATURALLY_NON_VEG), NATURALLY_NON_VEG := as.character(MISSCODE)]


  ## NON FORESTED ANTHROPOGENIC
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = NonFATable,
                            dtVar = "ANTH_NON", correspVar = "AVI",
                            newVar = "CASFRI", newName = "NON_FORESTED_ANTHRO", keepOld = TRUE)
  casfriDT[!is.na(ANTH_NON) & is.na(NON_FORESTED_ANTHRO), NON_FORESTED_ANTHRO := as.character(ERRCODE)]
  casfriDT[is.na(NON_FORESTED_ANTHRO), NON_FORESTED_ANTHRO := as.character(MISSCODE)]

  ## if NATURALLY_NON_VEG resulted in ERRCODE and NON_FORESTED_ANTHRO was not resolved, try matching
  ## the NAT_NON AVI code with NON_FORESTED_ANTHRO table
  temp <- invent2CASFRI(dt = casfriDT[NATURALLY_NON_VEG %in% ERRCODE], correspTab = NonFATable,
                        dtVar = "NAT_NON", correspVar = "AVI",
                        newVar = "CASFRI", newName = "NON_FORESTED_ANTHRO", keepOld = TRUE)
  casfriDT[NATURALLY_NON_VEG %in% as.character(ERRCODE) & is.na(NON_FORESTED_ANTHRO), NON_FORESTED_ANTHRO := temp[, NON_FORESTED_ANTHRO]]

  ## if NON_FORESTED_ANTHRO resulted in ERRCODE and NATURALLY_NON_VEG was not resolved, try matching
  ## the NON_FORESTED_ANTHRO AVI code with NATURALLY_NON_VEG table
  temp <- invent2CASFRI(dt = casfriDT[NON_FORESTED_ANTHRO %in% ERRCODE], correspTab = NatNVTable,
                        dtVar = "ANTH_NON", correspVar = "AVI",
                        newVar = "CASFRI", newName = "NATURALLY_NON_VEG", keepOld = TRUE)
  casfriDT[NON_FORESTED_ANTHRO %in% as.character(ERRCODE) & is.na(NATURALLY_NON_VEG), NATURALLY_NON_VEG := temp[, NATURALLY_NON_VEG]]


  ## NON FORESTED VEGETATION
  ## adjust values before conversion
  casfriDT[, NFL := NFLAdjustAB(nonForestLand = NFL, sp1 = SP1, layerID = LAYER_RANK),
           by = "P_ID"]
  casfriDT[toupper(ANTH_NON) %in% "AIL" & !is.na(NFL), NFL := ANTH_NON]

  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = NonFVTable,
                            dtVar = "NFL", correspVar = "AVI",
                            newVar = "CASFRI", newName = "NON_FORESTED_VEG", keepOld = TRUE)
  casfriDT[NFL %in% c("SO", "SC") &
             !is.na(HEIGHT) & (HEIGHT >= 0 & HEIGHT < 2),
           NON_FORESTED_VEG := "SL"]
  casfriDT[NFL %in% c("SO", "SC") &
             !is.na(HEIGHT) & (HEIGHT >= 2),
           NON_FORESTED_VEG := "ST"]
  casfriDT[NFL %in% c("SO", "SC") &
             (is.na(HEIGHT) | HEIGHT < 0),
           NON_FORESTED_VEG := as.character(ERRCODE)]

  casfriDT[!is.na(NFL) & is.na(NON_FORESTED_VEG), NON_FORESTED_VEG := as.character(ERRCODE)]
  casfriDT[is.na(NON_FORESTED_VEG), NON_FORESTED_VEG := as.character(MISSCODE)]


  ## DISTURBANCE CODES
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                            dtVar = "MOD1", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST1", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                            dtVar = "MOD2", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST2", keepOld = TRUE)

  casfriDT[!is.na(MOD1) & is.na(DIST1), DIST1 := as.character(ERRCODE)]
  casfriDT[!is.na(MOD2) & is.na(DIST2), DIST2 := as.character(ERRCODE)]
  casfriDT[is.na(MOD1), DIST1 := as.character(MISSCODE)]
  casfriDT[is.na(MOD2), DIST2 := as.character(MISSCODE)]

  if ("MOD3" %in% names(casfriDT)) {
    casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                              dtVar = "MOD3", correspVar = "AVI",
                              newVar = "CASFRI", newName = "DIST3", keepOld = TRUE)

    casfriDT[!is.na(MOD3) & is.na(DIST3), DIST3 := as.character(ERRCODE)]
    casfriDT[is.na(MOD3), DIST3 := as.character(MISSCODE)]
  } else {
    casfriDT[LAYER_RANK == 1, DIST3 := as.character(MISSCODE)]
    casfriDT[LAYER_RANK == 2, DIST3 := as.character(UNDEF)]   ## not sure why, but how it's coded in CASFRI perl script
  }


  ## DISTURBANCE YEARS
  casfriDT[, DIST1_YEAR := MOD1_YR]
  casfriDT[, DIST2_YEAR := MOD2_YR]
  casfriDT[MOD1_YR %in% 0, DIST1_YEAR :=  MISSCODE]
  casfriDT[MOD2_YR %in% 0, DIST2_YEAR :=  MISSCODE]
  casfriDT[MOD1_YR < 1880 & MOD1_YR > 2020, DIST1_YEAR :=  ERRCODE]
  casfriDT[MOD2_YR < 1880 & MOD2_YR > 2020, DIST1_YEAR :=  ERRCODE]

  if ("MOD3_YR" %in% names(casfriDT)) {
    casfriDT[, DIST3_YEAR := MOD3_YR]
    casfriDT[MOD3_YR %in% 0, DIST3_YEAR :=  MISSCODE]
    casfriDT[MOD3_YR < 1880 & MOD3_YR > 2020, DIST3_YEAR :=  ERRCODE]
  } else {
    casfriDT[LAYER_RANK == 1, DIST3_YEAR := as.character(MISSCODE)]
    casfriDT[LAYER_RANK == 2, DIST3_YEAR := as.character(UNDEF)]
  }

  ## DISTURBANCE EXTENSION UPPER AND LOWER
  ## UPPER
  casfriDT[, `:=` (MOD1_EXT = as.integer(MOD1_EXT),
                   MOD2_EXT = as.integer(MOD2_EXT))]
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                            dtVar = "MOD1_EXT", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST1_EXTENT_UPPER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                            dtVar = "MOD2_EXT", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST2_EXTENT_UPPER", keepOld = TRUE)

  ## LOWER
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                            dtVar = "MOD1_EXT", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST1_EXTENT_LOWER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                            dtVar = "MOD2_EXT", correspVar = "AVI",
                            newVar = "CASFRI", newName = "DIST2_EXTENT_LOWER", keepOld = TRUE)


  casfriDT[!is.na(MOD1_EXT) & is.na(DIST1_EXTENT_UPPER), DIST1_EXTENT_UPPER := ERRCODE]
  casfriDT[!is.na(MOD2_EXT) & is.na(DIST2_EXTENT_UPPER), DIST2_EXTENT_UPPER := ERRCODE]
  casfriDT[is.na(MOD1_EXT), DIST1_EXTENT_UPPER := MISSCODE]
  casfriDT[is.na(MOD2_EXT), DIST2_EXTENT_UPPER := MISSCODE]
  casfriDT[!is.na(MOD1_EXT) & is.na(DIST1_EXTENT_LOWER), DIST1_EXTENT_LOWER := ERRCODE]
  casfriDT[!is.na(MOD2_EXT) & is.na(DIST2_EXTENT_LOWER), DIST2_EXTENT_LOWER := ERRCODE]
  casfriDT[is.na(MOD1_EXT), DIST1_EXTENT_LOWER := MISSCODE]
  casfriDT[is.na(MOD2_EXT), DIST2_EXTENT_LOWER := MISSCODE]


  if ("MOD3_EXT" %in% names(casfriDT)) {
    casfriDT[, MOD3_EXT := as.integer(MOD3_EXT)]
    casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                              dtVar = "MOD3_EXT", correspVar = "AVI",
                              newVar = "CASFRI", newName = "DIST3_EXTENT_UPPER", keepOld = TRUE)
    casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                              dtVar = "MOD3_EXT", correspVar = "AVI",
                              newVar = "CASFRI", newName = "DIST3_EXTENT_LOWER", keepOld = TRUE)
    casfriDT[!is.na(MOD3_EXT) & is.na(DIST3_EXTENT_UPPER), DIST3_EXTENT_UPPER := as.character(ERRCODE)]
    casfriDT[is.na(MOD3_EXT), DIST3_EXTENT_UPPER := as.character(MISSCODE)]
    casfriDT[!is.na(MOD3_EXT) & is.na(DIST3_EXTENT_LOWER), DIST3_EXTENT_LOWER := as.character(ERRCODE)]
    casfriDT[is.na(MOD3_EXT), DIST3_EXTENT_LOWER := as.character(MISSCODE)]
  } else {
    casfriDT[LAYER_RANK == 1, DIST3_EXTENT_UPPER := as.character(MISSCODE)]
    casfriDT[LAYER_RANK == 1, DIST3_EXTENT_LOWER := as.character(MISSCODE)]
    casfriDT[LAYER_RANK == 2, DIST3_EXTENT_UPPER := as.character(UNDEF)]
    casfriDT[LAYER_RANK == 2, DIST3_EXTENT_LOWER := as.character(UNDEF)]
  }


  ## adjustments
  casfriDT[DIST1 %in% as.character(ERRCODE) &
             MOD1 %in% illegalDisturbances &
             (!DIST1_YEAR %in% c(MISSCODE, ERRCODE) | !is.na(DIST1_EXTENT_UPPER)),
           DIST1 := "OT"]
  casfriDT[DIST2 %in% as.character(ERRCODE) &
             MOD2 %in% illegalDisturbances &
             (!DIST2_YEAR %in% c(MISSCODE, ERRCODE) | !is.na(DIST2_EXTENT_UPPER)),
           DIST2 := "OT"]

  if ("MOD3_EXT" %in% names(casfriDT)) {
    casfriDT[DIST3 %in% as.character(ERRCODE) &
               MOD3 %in% illegalDisturbances &
               (!DIST3_YEAR %in% c(MISSCODE, ERRCODE) | !is.na(DIST3_EXTENT_UPPER)),
             DIST3 := "OT"]
  }

  casfriDT[toupper(ANTH_NON) %in% "AIL",
           `:=` (DIST1 = "CO", DIST1_EXTENT_UPPER = 100, DIST1_EXTENT_LOWER = 95)]

  ## PRODUCTIVE FOREST
  ## not present in the CAFRI manual, but present in perl scripts (both for SK an AB)
  ## in AB productive code is computed for canopy and understory layers alike
  casfriDT[, PRODUCTIVE_FOREST := "PF"]
  casfriDT[is.na(SP1) &
             (((!CROWN_CLOSURE_UPPER %in% MISSCODE | !CROWN_CLOSURE_LOWER %in% MISSCODE) &
                 !is.na(DENSITY)) | !HEIGHT_UPPER %in% MISSCODE | !HEIGHT_LOWER %in% MISSCODE),
           PRODUCTIVE_FOREST := "PP"]
  ## adjust species for PP
  spCols <- grep("SPEC[[:digit:]]*$", names(casfriDT), value = TRUE)
  casfriDT[PRODUCTIVE_FOREST %in% "PP",
           grep("SPEC[[:digit:]]*$", names(casfriDT), value = TRUE) :=
             as.list(rep(SPECIES_MISSCODE, length(spCols)))]

  spCols <- grep("SPEC.*_PER$", names(casfriDT), value = TRUE)
  casfriDT[PRODUCTIVE_FOREST %in% "PP",
           grep("SPEC.*_PER$", names(casfriDT), value = TRUE) :=
             as.list(rep(0, length(spCols)))]

  casfriDT[LAYER_RANK == 2 & DIST1 %in% "CO", PRODUCTIVE_FOREST := "PF"]

  ## WETLAND - here apply and test function
  casfriDT$WETLAND_CLASS <- as.character()
  casfriDT$WETLAND_VEG_MOD <- as.character()
  casfriDT$WETLAND_LAND_MOD <- as.character()
  casfriDT$WETLAND_LOCAL_MOD <- as.character()

  casfriDT[, c("WETLAND_CLASS", "WETLAND_VEG_MOD",
               "WETLAND_LAND_MOD", "WETLAND_LOCAL_MOD") :=
             wetlandCodesAB(SMR, DENSITY, NFL, NAT_NON, SP1, SP2,
                            SPEC1_PER, LAYER_RANK,
                            MISSCODE = MISSCODE, ERRCODE = ERRCODE),
           by = P_ID]

  ## add remaining variables that are not defined/missing for Saskatchewan
  set(casfriDT, NULL, setdiff(casfriVars, names(casfriDT)),
      value = UNDEF)

  ## remove old variables
  set(casfriDT, NULL, setdiff(names(casfriDT[, -"P_ID"]), casfriVars),
      value = NULL)

  ## CONVERT TO SF OBJECT
  casfriDT[order(P_ID)]
  if (any(!casfriDT$P_ID %in% inv$P_ID)) {
    stop(paste0("Polygon order doesn't match between CASFRI data.table and '", inv, "'"))
  } else
    invOut <- st_set_geometry(casfriDT, st_geometry(inv))

  return(invOut)
}

## Saskatchewan
SKToCASFRI <- function(inv, tablesDir, folder, dim) {
  message("Converting Saskatchewan inventory data to CASFRI standard...")
  ## CASFRI VARIABLES
  casfriVars <- c("CROWN_CLOSURE_LOWER", "CROWN_CLOSURE_UPPER", "DIST1",
                  "DIST1_EXTENT_LOWER", "DIST1_EXTENT_UPPER", "DIST1_YEAR",
                  "DIST2", "DIST2_EXTENT_LOWER", "DIST2_EXTENT_UPPER",
                  "DIST2_YEAR", "DIST3", "DIST3_EXTENT_LOWER",
                  "DIST3_EXTENT_UPPER", "DIST3_YEAR", "ECOSITE",
                  "HEIGHT_LOWER", "HEIGHT_UPPER", "LAYER", "LAYER_RANK",
                  "MAPSHEET_ID", "NATURALLY_NON_VEG", "NON_FORESTED_ANTHRO",
                  "NON_FORESTED_VEG", "NUMBER_OF_LAYERS", "ORIGIN_LOWER",
                  "ORIGIN_UPPER", "ORIGINAL_STAND_ID", "PHOTO_YEAR_MAX",
                  "PHOTO_YEAR_MIN", "POLYGON_AREA", "POLYGON_PERIMITER",
                  "PRODUCTIVE_FOREST", "SITE_CLASS", "SITE_INDEX", "SMR",
                  "SPEC1", "SPEC1_PER",
                  "SPEC2", "SPEC2_PER", "SPEC3",
                  "SPEC3_PER", "SPEC4", "SPEC4_PER",
                  "SPEC5", "SPEC5_PER", "SPEC6",
                  "SPEC6_PER", "SPEC7", "SPEC7_PER", "SPEC8",
                  "SPEC8_PER","SPEC9", "SPEC9_PER","SPEC10",
                  "SPEC10_PER", "STAND_STRUCTURE", "STAND_STRUCTURE_PER",
                  "STAND_STRUCTURE_RANGE", "UNPRODUCTIVE_FOREST", "WETLAND_CLASS",
                  "WETLAND_LAND_MOD", "WETLAND_LOCAL_MOD", "WETLAND_VEG_MOD")

  ## CASFRI ERROR/MISSING CODES
  # MISSCODE <- -1111L
  # ERRCODE <- -9999L
  # SPECIES_MISSCODE <- "XXXX MISS"
  # SPECIES_ERRCODE <- "XXXX ERRC"
  # UNDEF <- -8888L
  #
  MISSCODE <- NA
  ERRCODE <- NA
  SPECIES_MISSCODE <- "NA"
  SPECIES_ERRCODE <- "NA"
  UNDEF <- NA

  ## get all necessary conversion tables
  SMRTable <- read.xlsx(tablesDir, sheetName = "SMRTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)
  SStrucTable <- read.xlsx(tablesDir, sheetName = "StandStructTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)
  SppTable <- read.xlsx(tablesDir, sheetName = "SpeciesTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)
  distTable <- read.xlsx(tablesDir, sheetName = "DisturbanceTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)
  distExtUTable <- read.xlsx(tablesDir, sheetName = "DistExtUpperTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)
  distExtLTable <- read.xlsx(tablesDir, sheetName = "DistExtLowerTable", header = TRUE) %>%
    data.table(.) %>%
    .[, .(CASFRI, SFVI)] %>%
    na.omit(.)


  casfriDT <- st_set_geometry(inv, NULL) %>%
    data.table(.)

  ## MAPSHEET
  casfriDT[, MAPSHEET_ID := MAPSHEET_NUM]

  ## LAYER
  casfriDT[!is.na(LAYER),
           LAYER := as.character(as.numeric(factor(LAYER, levels = c("1", "2", "3", "S", "H"),
                                                   labels =  c("1", "2", "3", "S", "H")))),
           by = P_ID]
  casfriDT$LAYER <- as.numeric(casfriDT$LAYER)  ## can't coerce directly to numeric above.

  ## LAYER RANK
  ## 1, 2, 3 for crown layers, 4 for shrub and 5 for herbs
  casfriDT[, LAYER_RANK := LAYER]


  ## SOIL MOISTURE REGIME
  setnames(casfriDT, old = "SMR", new = "SMRsfvi") ## change name of original variable
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = SMRTable,
                            dtVar = "SMRsfvi", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "SMR", keepOld = TRUE)
  ## change NAs to CASFRI missing code
  casfriDT[is.na(SMRsfvi), SMR := as.character(MISSCODE)]
  ## replace NAs generated from unknown codes CASFRI by error code
  casfriDT[!is.na(SMRsfvi) & is.na(SMR), SMR := as.character(ERRCODE)]

  ## STAND STRUCTURE
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = SStrucTable,
                            dtVar = "LAYER_TYPE", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "STAND_STRUCTURE", keepOld = TRUE)
  ## replace unknown codes with errorcode
  casfriDT[!is.na(LAYER_TYPE) & !LAYER_TYPE %in% c("S", "C", "M"), STAND_STRUCTURE := as.character(ERRCODE)]
  ## change NAs to "S"
  casfriDT[is.na(LAYER_TYPE), STAND_STRUCTURE := "S"]

  ## NUMBER OF LAYERS
  ## in CASFRI num. layers = no. non-NA layer types per stand.
  ## for this dataset it is best to cound number of non-NA sp1 per stand
  casfriDT[, NUMBER_OF_LAYERS := sum(!is.na(LAYER_TYPE) & LAYER_TYPE != ERRCODE), by = P_ID]


  ## CANOPY CLOSURE UPPER & LOWER
  ## same as CC for Saskatchewan
  casfriDT[CROWN_CLOSURE >= 0 & CROWN_CLOSURE <= 100, CROWN_CLOSURE_UPPER := CROWN_CLOSURE]
  casfriDT[CROWN_CLOSURE >= 0 & CROWN_CLOSURE <= 100, CROWN_CLOSURE_LOWER := CROWN_CLOSURE]
  ## change NAs to missing code
  casfriDT[is.na(CROWN_CLOSURE), `:=` (CROWN_CLOSURE_UPPER = MISSCODE,
                                       CROWN_CLOSURE_LOWER = MISSCODE)]
  ## replace unknown values CASFRI by error code
  casfriDT[CROWN_CLOSURE < 0 | CROWN_CLOSURE > 100,`:=` (CROWN_CLOSURE_UPPER = ERRCODE,
                                                         CROWN_CLOSURE_LOWER = ERRCODE)]

  ## STAND HEIGHT UPPER & LOWER
  casfriDT[HEIGHT > 0 & HEIGHT <= 100, HEIGHT_UPPER:= HEIGHT + 0.5]
  casfriDT[HEIGHT >= 0.5 & HEIGHT <= 100, HEIGHT_LOWER := HEIGHT - 0.5]
  casfriDT[HEIGHT > 0 & HEIGHT < 0.5, HEIGHT_LOWER := HEIGHT]
  ## change NAs to missing code
  casfriDT[is.na(HEIGHT), `:=` (HEIGHT_UPPER = MISSCODE,
                                HEIGHT_LOWER = MISSCODE)]
  ## replace unknown values CASFRI by error code
  casfriDT[HEIGHT < 0 | HEIGHT > 100,`:=` (HEIGHT_UPPER = ERRCODE,
                                           HEIGHT_LOWER = ERRCODE)]
  ## special case where error code comes from missing layers
  casfriDT[LAYER %in% c(1,2,3) & (HEIGHT_LOWER %in% ERRCODE | HEIGHT_UPPER %in% ERRCODE)] %>%
    .[is.na(SP1) & HEIGHT %in% 0, `:=` (HEIGHT_UPPER = MISSCODE,
                                      HEIGHT_LOWER = MISSCODE)]

  ## SPECIES - CANOPY
  ## note that these fields will be empty for non canopy layers
  casfriDT[LAYER %in% c(1:3),
           `:=` (SPEC1 = spLatinName(SP1, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC2 = spLatinName(SP2, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC3 = spLatinName(SP3, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC4 = spLatinName(SP4, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC5 = spLatinName(SP5, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC6 = spLatinName(SP6, SppTable, MISSCODE = SPECIES_MISSCODE, ERRCODE = SPECIES_ERRCODE),
                 SPEC7 = as.character(UNDEF),
                 SPEC8 = as.character(UNDEF),
                 SPEC9 = as.character(UNDEF),
                 SPEC10 = as.character(UNDEF)),
           by = seq_len(nrow(casfriDT[LAYER %in% c(1:3)]))]

  ## SPECIES PERCENT - CANOPY
  ## note that these fields will be empty for non canopy layers
  casfriDT[LAYER %in% c(1:3),
           `:=` (SPEC1_PER = spPercent(SP1_COVER, province = "SK"),
                 SPEC2_PER = spPercent(SP2_COVER, province = "SK"),
                 SPEC3_PER = spPercent(SP3_COVER, province = "SK"),
                 SPEC4_PER = spPercent(SP4_COVER, province = "SK"),
                 SPEC5_PER = spPercent(SP5_COVER, province = "SK"),
                 SPEC6_PER = spPercent(SP6_COVER, province = "SK"),
                 SPEC7_PER = as.character(UNDEF),
                 SPEC8_PER = as.character(UNDEF),
                 SPEC9_PER = as.character(UNDEF),
                 SPEC10_PER = as.character(UNDEF)),
           by = seq_len(nrow(casfriDT[LAYER %in% c(1:3)]))]

  # casfriDT[sum(SP1_COVER, SP2_COVER, SP3_COVER, SP4_COVER, SP5_COVER) %in% 90,
  #          SPEC6_PER := 10]
  casfriDT[LAYER %in% c(1:3), c("SPEC1_PER", "SPEC2_PER", "SPEC3_PER",
                                "SPEC4_PER", "SPEC5_PER", "SPEC6_PER") :=
             spPercentAdjust(sppPer = c(SPEC1_PER,
                                        SPEC2_PER,
                                        SPEC3_PER,
                                        SPEC4_PER,
                                        SPEC5_PER,
                                        SPEC6_PER),
                             province = "SK"),
           by = seq_len(nrow(casfriDT[LAYER %in% c(1:3)]))]

  ## ORIGIN UPPER & LOWER
  casfriDT[,  `:=` (ORIGIN_UPPER = as.integer(originUpper(YOO, province= "SK", MISSCODE = MISSCODE, ERRCODE = ERRCODE)),
                    ORIGIN_LOWER = as.integer(originLower(YOO, province= "SK", MISSCODE = MISSCODE, ERRCODE = ERRCODE))),
           by = 1:nrow(casfriDT)]

  casfriDT[LAYER_RANK %in% c(1,2,3) &
             (ORIGIN_UPPER %in% ERRCODE | ORIGIN_LOWER %in% ERRCODE) &
             is.na(SP1) &
             (YOO %in% 0),
           `:=` (ORIGIN_UPPER = MISSCODE, ORIGIN_LOWER = MISSCODE)]

  ## NATURALLY_NON_VEG
  ## nonVegNatSK follows SK_conversion06.pm and the CASFRI manual, but because aquatic class has unknown classes,
  ## I'll use a modification of TYPE
  ## about non forested polygons
  casfriDT$NATURALLY_NON_VEG <-as.character()
  casfriDT[TYPE == "WAT", AQUATIC_CLASS2 := "WA"]
  casfriDT[LAYER %in% c(1,2,3),
           NATURALLY_NON_VEG := nonVegNatSK(nonForest = NVSL,
                                            aquatic = AQUATIC_CLASS2,
                                            ERRCODE = as.character(ERRCODE), MISSCODE = as.character(MISSCODE)),
           by = P_ID]

  ## NON FORESTED VEGETATION
  casfriDT$NON_FORESTED_VEG <- as.character()
  casfriDT[, NON_FORESTED_VEG := nonForestVegSK(sp1 = SP1,
                                                cover = CROWN_CLOSURE,
                                                moist = SMRsfvi,
                                                layerID = LAYER_RANK,
                                                ERRCODE = ERRCODE, MISSCODE = MISSCODE),
           by = P_ID]

  ## NON FORESTED ANTHROPOGENIC - LUC/NVSL empty in Dave's data.
  casfriDT$NON_FORESTED_ANTHRO <- as.character()
  casfriDT[, NON_FORESTED_ANTHRO := nonVegAnthSK(landUse = LUC,
                                                 nonForest = NVSL,
                                                 ERRCODE = ERRCODE, MISSCODE = MISSCODE),
           by = P_ID]


  ## if all got ERRCODE replace by MISSCODE
  casfriDT[NATURALLY_NON_VEG %in% as.character(ERRCODE) &
             NON_FORESTED_VEG %in% as.character(ERRCODE) &
             NON_FORESTED_ANTHRO %in% as.character(ERRCODE),
           c("NATURALLY_NON_VEG", "NON_FORESTED_VEG", "NON_FORESTED_ANTHRO") :=
             list(as.character(MISSCODE), as.character(MISSCODE), as.character(MISSCODE))]


  ## DISTURBANCE CODES
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                            dtVar = "DISTURBANCE_1", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST1", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                            dtVar = "DISTURBANCE_2", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST2", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distTable,
                            dtVar = "DISTURBANCE_3", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST3", keepOld = TRUE)
  casfriDT[is.na(DISTURBANCE_1), DIST1 := as.character(MISSCODE)]
  casfriDT[is.na(DISTURBANCE_2), DIST2 := as.character(MISSCODE)]
  casfriDT[is.na(DISTURBANCE_3), DIST3 := as.character(MISSCODE)]

  ## DISTURBANCE YEARS
  casfriDT[, DIST1_YEAR := YOD_1]
  casfriDT[, DIST2_YEAR := YOD_2]
  casfriDT[, DIST3_YEAR := YOD_3]
  casfriDT[YOD_1 %in% 0, DIST1_YEAR :=  MISSCODE]
  casfriDT[YOD_2 %in% 0, DIST2_YEAR :=  MISSCODE]
  casfriDT[YOD_3 %in% 0, DIST3_YEAR :=  MISSCODE]

  ## DISTURBANCE EXTENSION UPPER AND LOWER
  casfriDT[, `:=` (DISTURBANCE_EXTENT_1 = as.numeric(DISTURBANCE_EXTENT_1),
                   DISTURBANCE_EXTENT_2 = as.numeric(DISTURBANCE_EXTENT_2),
                   DISTURBANCE_EXTENT_3 = as.numeric(DISTURBANCE_EXTENT_3))]

  ## UPPER
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                            dtVar = "DISTURBANCE_EXTENT_1", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST1_EXTENT_UPPER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                            dtVar = "DISTURBANCE_EXTENT_2", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST2_EXTENT_UPPER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtUTable,
                            dtVar = "DISTURBANCE_EXTENT_3", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST3_EXTENT_UPPER", keepOld = TRUE)

  ## LOWER
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                            dtVar = "DISTURBANCE_EXTENT_1", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST1_EXTENT_LOWER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                            dtVar = "DISTURBANCE_EXTENT_2", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST2_EXTENT_LOWER", keepOld = TRUE)
  casfriDT <- invent2CASFRI(dt = casfriDT, correspTab = distExtLTable,
                            dtVar = "DISTURBANCE_EXTENT_3", correspVar = "SFVI",
                            newVar = "CASFRI", newName = "DIST3_EXTENT_LOWER", keepOld = TRUE)

  ## deal with missing values
  casfriDT[is.na(DISTURBANCE_EXTENT_1), DIST1_EXTENT_UPPER := MISSCODE]
  casfriDT[is.na(DISTURBANCE_EXTENT_2), DIST2_EXTENT_UPPER := MISSCODE]
  casfriDT[is.na(DISTURBANCE_EXTENT_3), DIST3_EXTENT_UPPER := MISSCODE]
  casfriDT[is.na(DISTURBANCE_EXTENT_1), DIST1_EXTENT_LOWER := MISSCODE]
  casfriDT[is.na(DISTURBANCE_EXTENT_2), DIST2_EXTENT_LOWER := MISSCODE]
  casfriDT[is.na(DISTURBANCE_EXTENT_3), DIST3_EXTENT_LOWER := MISSCODE]

  ## NAs generated after conversion get ERRCODE
  casfriDT[is.na(DIST1_EXTENT_UPPER), DIST1_EXTENT_UPPER := ERRCODE]
  casfriDT[is.na(DIST1_EXTENT_UPPER), DIST2_EXTENT_UPPER := ERRCODE]
  casfriDT[is.na(DIST1_EXTENT_UPPER), DIST3_EXTENT_UPPER := ERRCODE]
  casfriDT[is.na(DIST1_EXTENT_LOWER), DIST1_EXTENT_LOWER := ERRCODE]
  casfriDT[is.na(DIST1_EXTENT_LOWER), DIST2_EXTENT_LOWER := ERRCODE]
  casfriDT[is.na(DIST1_EXTENT_LOWER), DIST3_EXTENT_LOWER := ERRCODE]

  ## PRODUCTIVE FOREST
  ## not present in the CASFRI manual, but present in perl scripts (both for SK an AB)
  ## this doesn't make sens because it ends up classifying water polygons as productive forest
  casfriDT[, PRODUCTIVE_FOREST := "PF"]
  casfriDT[LAYER_RANK == 1 & is.na(SP1) &
             (CROWN_CLOSURE_LOWER > 0 | HEIGHT_LOWER > 0) &
             any(!DISTURBANCE_1 %in% "CO" & !DISTURBANCE_2 %in% "CO" & !DISTURBANCE_3 %in% "CO"),
           PRODUCTIVE_FOREST := "PP"]


  ## WETLAND
  ## wetlandCodesSK follows the approach in SK_conversion06MISTIK.pm
  ## wetlandCodesSK2 is an adaptation because Dave's data is empty on NVSL
  casfriDT$WETLAND_CLASS <- as.character()
  casfriDT$WETLAND_VEG_MOD <- as.character()
  casfriDT$WETLAND_LAND_MOD <- as.character()
  casfriDT$WETLAND_LOCAL_MOD <- as.character()

  # casfriDT[, c("WETLAND_CLASS", "WETLAND_VEG_MOD",
  #              "WETLAND_LAND_MOD", "WETLAND_LOCAL_MOD") :=
  #            wetlandCodesSK(SMRsfvi, CROWN_CLOSURE, HEIGHT, NVSL, SP1, SP2,
  #                           SPEC1_PER, LAYER_RANK,
  #                           MISSCODE = MISSCODE, ERRCODE = ERRCODE),
  #          by = P_ID]

  casfriDT[, c("WETLAND_CLASS", "WETLAND_VEG_MOD",
               "WETLAND_LAND_MOD", "WETLAND_LOCAL_MOD") :=
             wetlandCodesSK2(SMRsfvi, CROWN_CLOSURE, HEIGHT, TYPE, SP1, SP2,
                             SPEC1_PER, LAYER_RANK,
                             MISSCODE = MISSCODE, ERRCODE = ERRCODE),
           by = P_ID]

  ## add remaining variables that are not defined/missing for Saskatchewan
  set(casfriDT, NULL, setdiff(casfriVars, names(casfriDT)),
      value = UNDEF)

  ## remove old variables
  set(casfriDT, NULL, setdiff(names(casfriDT[, -"P_ID"]), casfriVars),
      value = NULL)

  ## CONVERT TO SF OBJECT
  casfriDT[order(P_ID)]
  if (any(!casfriDT$P_ID %in% inv$P_ID)) {
    stop(paste0("Polygon order doesn't match between CASFRI data.table and '", inv, "'"))
  } else
    invOut <- st_set_geometry(casfriDT, st_geometry(inv))

  return(invOut)
}

## MELT PRE-FIRE DATA FUNCTIONS -------
## these functions transform the Alberta/Saskatchewan pre-fire data from
## an extended table format to a molten format before the data can be converted to
## CASFRI standards
## inv is a shapefile containing inventory polygons
## allVars is a vector of final variables that should be present in the molten data (only for Alberta).
## folder is the directory path where the molten shapefiles will be saved
## dim is used for faster caching - omit "inv"
## Alberta
meltPreFireABInv <- function(inv, invName, allVars, folder, dim) {
  tmpDT <- st_set_geometry(inv, NULL) %>%
    data.table(.)

  ## IF SOME VARIABLES ARE MISSING ADD THEM AS NAs
  if (any(!allVars %in% names(tmpDT))) {
    addCols <- allVars[!allVars %in% names(tmpDT)]
    tmpDT[, c(addCols) := list("NA", "NA")]
  }

  ## melt various fields
  moistDT <- melt(tmpDT[, .(P_ID, MOIST_REG, UMOIST_REG)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "MOIST_REG") %>%
    .[, LAYER := ifelse(LAYER == "MOIST_REG", "1", "2")]

  densityDT <- melt(tmpDT[, .(P_ID, DENSITY, UDENSITY)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "DENSITY") %>%
    .[, LAYER := ifelse(LAYER == "DENSITY", "1", "2")]

  heightDT <- melt(tmpDT[, .(P_ID, HEIGHT, UHEIGHT)], id.vars = "P_ID",
                   variable.name = "LAYER", value.name = "HEIGHT") %>%
    .[, LAYER := ifelse(LAYER == "HEIGHT", "1", "2")]

  sp1DT <- melt(tmpDT[, .(P_ID, SP1, USP1)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "SP1") %>%
    .[, LAYER := ifelse(LAYER == "SP1", "1", "2")]

  sp2DT <- melt(tmpDT[, .(P_ID, SP2, USP2)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "SP2") %>%
    .[, LAYER := ifelse(LAYER == "SP2", "1", "2")]

  sp3DT <- melt(tmpDT[, .(P_ID, SP3, USP3)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "SP3") %>%
    .[, LAYER := ifelse(LAYER == "SP3", "1", "2")]

  sp4DT <- melt(tmpDT[, .(P_ID, SP4, USP4)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "SP4") %>%
    .[, LAYER := ifelse(LAYER == "SP4", "1", "2")]

  sp5DT <- melt(tmpDT[, .(P_ID, SP5, USP5)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "SP5") %>%
    .[, LAYER := ifelse(LAYER == "SP5", "1", "2")]

  sp1_PERDT <- melt(tmpDT[, .(P_ID, SP1_PER, USP1_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "SP1_PER") %>%
    .[, LAYER := ifelse(LAYER == "SP1_PER", "1", "2")]

  sp2_PERDT <- melt(tmpDT[, .(P_ID, SP2_PER, USP2_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "SP2_PER") %>%
    .[, LAYER := ifelse(LAYER == "SP2_PER", "1", "2")]

  sp3_PERDT <- melt(tmpDT[, .(P_ID, SP3_PER, USP3_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "SP3_PER") %>%
    .[, LAYER := ifelse(LAYER == "SP3_PER", "1", "2")]

  sp4_PERDT <- melt(tmpDT[, .(P_ID, SP4_PER, USP4_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "SP4_PER") %>%
    .[, LAYER := ifelse(LAYER == "SP4_PER", "1", "2")]

  sp5_PERDT <- melt(tmpDT[, .(P_ID, SP5_PER, USP5_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "SP5_PER") %>%
    .[, LAYER := ifelse(LAYER == "SP5_PER", "1", "2")]

  strucDT <- melt(tmpDT[, .(P_ID, STRUC, USTRUC)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "STRUC") %>%
    .[, LAYER := ifelse(LAYER == "STRUC", "1", "2")]

  strucVDT <- melt(tmpDT[, .(P_ID, STRUC_VAL, USTRUC_VAL)], id.vars = "P_ID",
                   variable.name = "LAYER", value.name = "STRUC_VAL") %>%
    .[, LAYER := ifelse(LAYER == "STRUC_VAL", "1", "2")]

  originDT <- melt(tmpDT[, .(P_ID, ORIGIN, UORIGIN)], id.vars = "P_ID",
                   variable.name = "LAYER", value.name = "ORIGIN") %>%
    .[, LAYER := ifelse(LAYER == "ORIGIN", "1", "2")]

  tprDT <- melt(tmpDT[, .(P_ID, TPR, UTPR)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "TPR") %>%
    .[, LAYER := ifelse(LAYER == "TPR", "1", "2")]

  initialsDT <- melt(tmpDT[, .(P_ID, INITIALS, UINITIALS)], id.vars = "P_ID",
                     variable.name = "LAYER", value.name = "INITIALS") %>%
    .[, LAYER := ifelse(LAYER == "INITIALS", "1", "2")]

  nflDT <- melt(tmpDT[, .(P_ID, NFL, UNFL)], id.vars = "P_ID",
                variable.name = "LAYER", value.name = "NFL") %>%
    .[, LAYER := ifelse(LAYER == "NFL", "1", "2")]

  nfl_PERDT <- melt(tmpDT[, .(P_ID, NFL_PER, UNFL_PER)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "NFL_PER") %>%
    .[, LAYER := ifelse(LAYER == "NFL_PER", "1", "2")]

  natnonDT <- melt(tmpDT[, .(P_ID, NAT_NON, UNAT_NON)], id.vars = "P_ID",
                   variable.name = "LAYER", value.name = "NAT_NON") %>%
    .[, LAYER := ifelse(LAYER == "NAT_NON", "1", "2")]

  anthvegDT <- melt(tmpDT[, .(P_ID, ANTH_VEG, UANTH_VEG)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "ANTH_VEG") %>%
    .[, LAYER := ifelse(LAYER == "ANTH_VEG", "1", "2")]

  anthnonDT <- melt(tmpDT[, .(P_ID, ANTH_NON, UANTH_NON)], id.vars = "P_ID",
                    variable.name = "LAYER", value.name = "ANTH_NON") %>%
    .[, LAYER := ifelse(LAYER == "ANTH_NON", "1", "2")]

  mod1DT <- melt(tmpDT[, .(P_ID, MOD1, UMOD1)], id.vars = "P_ID",
                 variable.name = "LAYER", value.name = "MOD1") %>%
    .[, LAYER := ifelse(LAYER == "MOD1", "1", "2")]

  mod1EDT <- melt(tmpDT[, .(P_ID, MOD1_EXT, UMOD1_EXT)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "MOD1_EXT") %>%
    .[, LAYER := ifelse(LAYER == "MOD1_EXT", "1", "2")]

  mod1YDT <- melt(tmpDT[, .(P_ID, MOD1_YR, UMOD1_YR)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "MOD1_YR") %>%
    .[, LAYER := ifelse(LAYER == "MOD1_YR", "1", "2")]

  mod2DT <- melt(tmpDT[, .(P_ID, MOD2, UMOD2)], id.vars = "P_ID",
                 variable.name = "LAYER", value.name = "MOD2") %>%
    .[, LAYER := ifelse(LAYER == "MOD2", "1", "2")]

  mod2EDT <- melt(tmpDT[, .(P_ID, MOD2_EXT, UMOD2_EXT)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "MOD2_EXT") %>%
    .[, LAYER := ifelse(LAYER == "MOD2_EXT", "1", "2")]

  mod2YDT <- melt(tmpDT[, .(P_ID, MOD2_YR, UMOD2_YR)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "MOD2_YR") %>%
    .[, LAYER := ifelse(LAYER == "MOD2_YR", "1", "2")]

  tmpDT[, UDATA := "NA"]
  dataDT <- melt(tmpDT[, .(P_ID, DATA, UDATA)], id.vars = "P_ID",
                 variable.name = "LAYER", value.name = "DATA") %>%
    .[, LAYER := ifelse(LAYER == "DATA", "1", "2")]

  tmpDT[, UDATA_YR := "NA"]
  dataYDT <- melt(tmpDT[, .(P_ID, DATA_YR, UDATA_YR)], id.vars = "P_ID",
                  variable.name = "LAYER", value.name = "DATA_YR") %>%
    .[, LAYER := ifelse(LAYER == "DATA_YR", "1", "2")]

  ## make a list of all data.tables to join
  dtToJoin <- grep("^layer.*DT$|^density.*DT$|^anth.*DT$|^data.*DT$|^sp.*DT$|^height.*DT$|^origin.*DT$|^initials.*DT$|^mod.*DT$|^moist.*DT$|^nat.*DT$|^nfl.*DT$|^tpr.*DT$|^struc.*DT$",
                   ls(), value = TRUE)
  dtToJoin <- lapply(dtToJoin, get, envir = environment())

  ## join data.tables by P_ID and layer
  jointDT <- Reduce(function(x,y) merge(x, y, by = c("P_ID", "LAYER"), all = TRUE), dtToJoin)

  ## add other variables that dP_ID not need melting
  cols <- c(grep(paste0(names(jointDT), collapse = "|"),
                 names(tmpDT), value = TRUE, invert = TRUE),
            "P_ID")

  jointDT <- jointDT[tmpDT[, ..cols], on = "P_ID"]

  ## make new SF object with molten data
  invOutName <- paste0(invName, "Melt")

  invOut <- jointDT[data.table(P_ID = tmpDT$P_ID, geometry = inv[, 'geometry', drop = TRUE]), on = 'P_ID'] %>%
    as.data.frame(.) %>%
    st_sf(.)
  st_write(invOut, dsn = file.path(folder, paste0(invOutName, ".shp")), delete_layer = TRUE)

  return(invOut)
}

## Saskatechewan
meltPreFireSKInv <- function(inv,  invName, folder, dim) {
  tmpDT <- st_set_geometry(inv, NULL) %>%
    data.table(.)

  ## melt various fields
  layerDT <- melt(tmpDT[, .(P_ID, LAYER1_, LAYER2_, LAYER3_)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "LAYER") %>%
    .[, LAYER_NUM := sub("LAYER(.*?)_", "\\1", sub("ERBS_C|HRUB_C", "", LAYER_NUM))]

  canopyDT <- melt(tmpDT[, .(P_ID, CANOPY1, CANOPY2, CANOPY3)], id.vars = "P_ID",
                   variable.name = "LAYER_NUM", value.name = "CANOPY") %>%
    .[, LAYER_NUM := sub("CANOPY", "", LAYER_NUM)]

  crownDT <- melt(tmpDT[, .(P_ID, CROWN1_, CROWN2_, CROWN3_, HERBS_C, SHRUB_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "CROWN") %>%
    .[, LAYER_NUM := sub("CROWN(.*?)_", "\\1", sub("ERBS_C|HRUB_C", "", LAYER_NUM))]

  heightDT <- melt(tmpDT[, .(P_ID, HEIGHT1, HEIGHT2, HEIGHT3)], id.vars = "P_ID",
                   variable.name = "LAYER_NUM", value.name = "HEIGHT") %>%
    .[, LAYER_NUM := sub("HEIGHT", "", LAYER_NUM)]

  height_DT <- melt(tmpDT[, .(P_ID, HEIGHT1_, HEIGHT2_, HEIGHT3_)], id.vars = "P_ID",
                    variable.name = "LAYER_NUM", value.name = "HEIGHT_") %>%
    .[, LAYER_NUM := sub("HEIGHT(.*?)_", "\\1", LAYER_NUM)]

  sp1DT <- melt(tmpDT[, .(P_ID, HERBS1, SHRUB1, SP1_1, SP1_2, SP1_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP1") %>%
    .[, LAYER_NUM := sub("ERBS1|HRUB1|SP1_", "", LAYER_NUM)]

  sp2DT <- melt(tmpDT[, .(P_ID, HERBS2, SHRUB2, SP2_1, SP2_2, SP2_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP2") %>%
    .[, LAYER_NUM := sub("ERBS2|HRUB2|SP2_", "", LAYER_NUM)]

  sp3DT <- melt(tmpDT[, .(P_ID, HERBS3, SHRUB3, SP3_1, SP3_2, SP3_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP3") %>%
    .[, LAYER_NUM := sub("ERBS3|HRUB3|SP3_", "", LAYER_NUM)]

  sp4DT <- melt(tmpDT[, .(P_ID, HERBS4, SHRUB4, SP4_1, SP4_2, SP4_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP4") %>%
    .[, LAYER_NUM := sub("ERBS4|HRUB4|SP4_", "", LAYER_NUM)]

  sp5DT <- melt(tmpDT[, .(P_ID, HERBS5, SHRUB5, SP5_1, SP5_2, SP5_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP5") %>%
    .[, LAYER_NUM := sub("ERBS5|HRUB5|SP5_", "", LAYER_NUM)]

  sp6DT <- melt(tmpDT[, .(P_ID, HERBS6, SHRUB6, SP6_1, SP6_2, SP6_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP6") %>%
    .[, LAYER_NUM := sub("ERBS6|HRUB6|SP6_", "", LAYER_NUM)]

  sp7DT <- melt(tmpDT[, .(P_ID, HERBS7, SHRUB7)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP7") %>%
    .[, LAYER_NUM := sub("ERBS7|HRUB7", "", LAYER_NUM)]

  sp8DT <- melt(tmpDT[, .(P_ID, HERBS8, SHRUB8)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP8") %>%
    .[, LAYER_NUM := sub("ERBS8|HRUB8", "", LAYER_NUM)]

  sp9DT <- melt(tmpDT[, .(P_ID, HERBS9, SHRUB9)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "SP9") %>%
    .[, LAYER_NUM := sub("ERBS9|HRUB9", "", LAYER_NUM)]

  sp10DT <- melt(tmpDT[, .(P_ID, HERBS10, SHRUB10)], id.vars = "P_ID",
                 variable.name = "LAYER_NUM", value.name = "SP10") %>%
    .[, LAYER_NUM := sub("ERBS10|HRUB10", "", LAYER_NUM)]

  sp1_CDT <- melt(tmpDT[, .(P_ID, HERBS1_, SHRUB1_, SP1_1_C, SP1_2_C, SP1_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP1_C") %>%
    .[, LAYER_NUM := sub("SP1_(.*?)_C", "\\1", sub("ERBS1_|HRUB1_", "", LAYER_NUM))]

  sp2_CDT <- melt(tmpDT[, .(P_ID, HERBS2_, SHRUB2_, SP2_1_C, SP2_2_C, SP2_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP2_C") %>%
    .[, LAYER_NUM := sub("SP2_(.*?)_C", "\\1", sub("ERBS2_|HRUB2_", "", LAYER_NUM))]

  sp3_CDT <- melt(tmpDT[, .(P_ID, HERBS3_, SHRUB3_, SP3_1_C, SP3_2_C, SP3_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP3_C") %>%
    .[, LAYER_NUM := sub("SP3_(.*?)_C", "\\1", sub("ERBS3_|HRUB3_", "", LAYER_NUM))]

  sp4_CDT <- melt(tmpDT[, .(P_ID, HERBS4_, SHRUB4_, SP4_1_C, SP4_2_C, SP4_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP4_C") %>%
    .[, LAYER_NUM := sub("SP4_(.*?)_C", "\\1", sub("ERBS4_|HRUB4_", "", LAYER_NUM))]

  sp5_CDT <- melt(tmpDT[, .(P_ID, HERBS5_, SHRUB5_, SP5_1_C, SP5_2_C, SP5_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP5_C") %>%
    .[, LAYER_NUM := sub("SP5_(.*?)_C", "\\1", sub("ERBS5_|HRUB5_", "", LAYER_NUM))]

  sp6_cDT <- melt(tmpDT[, .(P_ID, HERBS6_, SHRUB6_, SP6_1_C, SP6_2_C, SP6_3_C)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP6_C") %>%
    .[, LAYER_NUM := sub("SP6_(.*?)_C", "\\1", sub("ERBS6_|HRUB6_", "", LAYER_NUM))]

  sp7_CDT <- melt(tmpDT[, .(P_ID, HERBS7_, SHRUB7_)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP7_C") %>%
    .[, LAYER_NUM := sub("ERBS7_|HRUB7_", "", LAYER_NUM)]

  sp8_CDT <- melt(tmpDT[, .(P_ID, HERBS8_, SHRUB8_)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP8_C") %>%
    .[, LAYER_NUM := sub("ERBS8_|HRUB8_", "", LAYER_NUM)]

  sp9_CDT <- melt(tmpDT[, .(P_ID, HERBS9_, SHRUB9_)], id.vars = "P_ID",
                  variable.name = "LAYER_NUM", value.name = "SP9_C") %>%
    .[, LAYER_NUM := sub("ERBS9_|HRUB9_", "", LAYER_NUM)]

  sp10_CDT <- melt(tmpDT[, .(P_ID, HERBS10_, SHRUB10_)], id.vars = "P_ID",
                   variable.name = "LAYER_NUM", value.name = "SP10_C") %>%
    .[, LAYER_NUM := sub("ERBS10_|HRUB10_", "", LAYER_NUM)]

  yooDT <- melt(tmpDT[, .(P_ID, YOO_1, YOO_2, YOO_3)], id.vars = "P_ID",
                variable.name = "LAYER_NUM", value.name = "YOO") %>%
    .[, LAYER_NUM := sub("YOO_", "", LAYER_NUM)]

  yoo_tyDT <- melt(tmpDT[, .(P_ID, YOO1_TY, YOO2_TY, YOO3_TY)], id.vars = "P_ID",
                   variable.name = "LAYER_NUM", value.name = "YOO_TY") %>%
    .[, LAYER_NUM := sub("YOO(.*?)_TY", "\\1", LAYER_NUM)]

  ## make a list of all data.tables to join
  dtToJoin <- grep("^layer.*DT$|^canopy.*DT$|^crown.*DT$|^height.*DT$|^sp.*DT$|^yoo.*DT$", ls(), value = TRUE)
  dtToJoin <- lapply(dtToJoin, get, envir = environment())

  ## join data.tables by P_ID and layer
  jointDT <- Reduce(function(x,y) merge(x, y, by = c("P_ID", "LAYER_NUM"), all = TRUE), dtToJoin)

  ## add other variables that dP_ID not need melting
  cols <- c(grep(paste0(paste0("^", c(names(jointDT),
                                      "HERB", "SHRUB"), ".*"), collapse = "|"),
                 names(tmpDT), value = TRUE, invert = TRUE),
            "P_ID")

  jointDT <- jointDT[tmpDT[, ..cols], on = "P_ID"]

  ## make new SF object with molten data
  invOutName <- paste0(invName, "Melt")

  invOut <- jointDT[data.table(P_ID = tmpDT$P_ID, geometry = inv[, 'geometry', drop = TRUE]), on = 'P_ID'] %>%
    as.data.frame(.) %>%
    st_sf(.)

  st_write(invOut, dsn = file.path(folder, paste0(invOutName, ".shp")), delete_layer = TRUE)

  return(invOut)
}


## DATA MISMATCH TROUBLESHOOTING FUNCTIONS -------
## these functions check for mismatches and correct data accordingly
## For now, functions were targeted for Saskatchewan data mismatches between the
## TYPE, CSG and PFT (final columns names, not original) colums

## add missing water info, if none other is available
addWaterInfo <- function(AQUATIC_CLASS, TYPE, CSG, PFT) {
  if (AQUATIC_CLASS > 0 & is.na(TYPE)) {
    if (is.na(CSG) & is.na(PFT))
      TYPE <- "WAT" else
        if (!CSG %in% "WAT" & !PFT %in% "WAT")
          warning("AQUATIC_CLASS is 'WAT', but CSG/PFT are not! No changes made")
  }
  TYPE
}

## correct mismatches between TYPE, CSG and PFT
correctCSGPFTTYPE <- function(LAYER, TYPE, CSG, PFT, SP1_COVER, SMR) {
  ## checks
  if (length(LAYER) > 5 |
      length(TYPE) > 5|
      length(CSG) > 5 |
      length(PFT) > 5|
      length(SP1_COVER) > 5 |
      length(SMR) > 5)
    stop("Saskatchewan: there should only be 5 different veg. layers")

  if(length(unique(TYPE)) > 1)
    stop("Saskatchewan: there should only be one TYPE value per stand")
  if(length(unique(CSG)) > 1)
    stop("Saskatchewan: there should only be one CSG value per stand")
  if(length(unique(PFT)) > 1)
    stop("Saskatchewan: there should only be one PFT value per stand")
  if(length(unique(SMR)) > 1)
    stop("Saskatchewan: there should only be one PFT value per stand")

  TYPE <- toupper(as.character(TYPE[LAYER == "1"]))
  CSG <- toupper(as.character(CSG[LAYER == "1"]))
  PFT <- toupper(as.character(PFT[LAYER == "1"]))
  SMR <- toupper(as.character(SMR[LAYER == "1"]))

  if (is.na(TYPE) | TYPE %in% "NFA") {
    if (CSG %in% c("WAT", "FOR", "TMS", "BSH", "OMS", "GRS", "RCK", "ALA", "UCL", "OTH") &
        PFT %in% c("WAT", "FOR", "TMS", "BSH", "OMS", "GRS", "RCK", "ALA", "UCL", "OTH")) {
      if (CSG %in% PFT) {
        TYPE <- CSG
        CSG <- "NA"
        PFT <- "NA"
      } else
        warning("CSG and PFT have a TYPE class, but differ. TYPE/CSG/PFT were not changed")
    } else {
      if (SP1_COVER[LAYER == "1"] > 0) {
        TYPE <- "FOR"
      } else if (sum(SP1_COVER) > 0) {
        if (SP1_COVER[LAYER == "S"] < SP1_COVER[LAYER == "H"]) {
          TYPE <- if (SMR %in% c("MW", "W", "VW")) "OMS" else "GRS"
        } else {
          TYPE <- if (SMR %in% c("MW", "W", "VW")) "BSH" else "OTH"
        }
      } else warning("Can't figure out TYPE class - no vegetation info")
    }
  }

  list(as.character(rep(TYPE, length(LAYER))),
       as.character(rep(CSG, length(LAYER))),
       as.character(rep(PFT, length(LAYER))))
}
