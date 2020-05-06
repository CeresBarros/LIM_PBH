#' Classification of stand structure into Cameron Naficy's vegetation types.
#' Uses a set of rules based on species relative biomass in a stand to
#' classify it  into one of 12 vegetation types:
#' "Oak", "PJ", "purePIPO", "DMCPIPO", "dryPSME", "PSME", "DMCPSME", "PICO",
#' "PIEN", "Broadleaf", "Mixedwood", "MMC".
#' This function deals with the eventuality of multiple matches between Cameron
#' species codes and simulated species codes
#'
#' @param Cameron character string of species in the stand according to Cameron's speces codes
#' @param speciesCode character string of species codes  that correspond to Cameron species -
#'    Note that because these are simualted species, one can correspond to several of Cameron's species
#'    e.g. Pinus sp can be PICO and PIFL
#' @param relB numeric vector of relative biomass for each species in each stand.
#' @param pureCutoff threshold of relative biomass above which a stand is considered pure. Defaults to
#'    0.8.
#' @param drySp character spring of species characteristic of dry sites. Defaults to
#'    \code{c("PSME", "PIPO", "PIFL", "JUSC", "QUGA")}
#' @param moistSp character spring of species characteristic of moist sites. Defaults to
#'    \code{c("ABLA", "BEPA", "PIEN", "PIGL", "PIMO", "POBA", "THPL")}
#'
#' @export

convertToCNVegType <- function (Cameron, speciesCode, relB, pureCutoff = 0.8,
                                drySp = c("PSME", "PIPO", "PIFL", "JUSC", "QUGA"),
                                moistSp = c("ABLA", "BEPA", "PIEN", "PIGL", "PIMO", "POBA", "THPL")) {
  tempDT <- data.table(Cameron, speciesCode, relB)
  ## Oak woodlands are dominated by oaks with no more dominant tree stature species
  oak <- all(.sumRelBs("QUGA", tempDT) >= pureCutoff,
             .sumRelBs(c('PIPO', 'PSME', 'PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS'), tempDT) < 0.05)

  ## P-J woodlands are dominated by Pinyon juniper trees with no more dominant tree stature species
  PJ <- all(.sumRelBs(c('PIED', 'PIMO2','JUSC', 'JUOC', 'JUOS', 'QUGA'), tempDT) >= pureCutoff,
            .sumRelBs(c('PIPO', 'PSME'), tempDT) < 0.05)

  ## Pure PIPO if PIPO is heavily dominant and accompanied by small amount of other species
  purePIPO <- all(.sumRelBs("PIPO", tempDT) >= pureCutoff,
                  .sumRelBs(c('PSME', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'), tempDT) < 0.30)

  ## DMC if PIPO present at >= 10% but less than 70% and other species are all dry site species
  DMCPIPO <- all(.sumRelBs("PIPO", tempDT) >= 0.10,
                 .sumRelBs("PIPO", tempDT) < pureCutoff,
                 .sumRelBs(c('PSME', 'PIPO', 'PIFL', 'PIED', 'PIMO2', 'JUSC', 'JUOC', 'JUOS', 'QUGA'), tempDT) >= 0.50)

  ## If PSME is dominant and dry site species are present
  dryPSME <- all(.sumRelBs("PSME", tempDT) >= pureCutoff,
                 .sumRelBs(moistSp, tempDT) < 0.10,
                 .sumRelBs(c('JUSC', 'JUOC', 'JUOS', 'PIFL', 'PIED', 'PIMO2', 'QUGA'), tempDT) > 0.05)
  ## If PSME is dominant and dry site species are absent
  PSME <- all(.sumRelBs("PSME", tempDT) >= pureCutoff,
              .sumRelBs(moistSp, tempDT) < 0.10,
              .sumRelBs(c('JUSC', 'JUOC', 'JUOS', 'PIFL', 'PIED', 'PIMO2', 'QUGA'), tempDT) <= 0.05)

  ## DMC if ponderosa pine not present, dry site species are dominant but may be
  ## micov.BAed with some other species (e.g. POTR, LAOC, PICO), and few moist site species are present in significant numbers
  DMCPSME <- all(.sumRelBs(drySp, tempDT) >= 0.50,
                 .sumRelBs(moistSp, tempDT) < 0.10)

  ## PICO if PICO dominates stand
  PICO <- .sumRelBs("PICO", tempDT) >= 0.5

  ## lowland PIEN if dominated by spruce
  PIEN <- .sumRelBs(c('PIEN', 'PIEN/PIGL', 'PIGL', 'ABLA'), tempDT) > 0.50

  ## broadleaf and mixedwood
  broadleaf <- .sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), tempDT) >= pureCutoff
  mixedwood <- all(.sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), tempDT) < pureCutoff,
                   .sumRelBs(c('POTR', 'POTR5', 'POBA', 'BEPA'), tempDT) >= 0.25)

  ## the final veg. type needs to be accessed sequentially.
  if (oak) return("Oak") else
    if (PJ) return("PJ") else
      if (purePIPO) return("purePIPO") else
        if (DMCPIPO) return("DMCPIPO") else
          if (dryPSME) return("dryPSME") else
            if (PSME) return("PSME") else
              if (DMCPSME) return("DMCPSME") else
                if (PICO) return("PICO") else
                  if(PIEN) return("PIEN") else
                    if (broadleaf) return("Broadleaf") else
                      if (mixedwood) ("Mixedwood") else
                        return("MMC")
}

#' internal function that sums relative biomasses for species matching a character string,
#' but that can be appear duplicated in another species coding column.
#' @param sppToMatch character string of species to match against for summing B.
#' @param DT data.table with columns 'Cameron', 'speciesCode', 'relB'.
.sumRelBs <- function(sppToMatch, DT) {
  DT[Cameron %in% sppToMatch] %>%
    .[, .(speciesCode, relB)] %>%
    unique(.) %>%
    .$relB %>%
    sum(.)
}
