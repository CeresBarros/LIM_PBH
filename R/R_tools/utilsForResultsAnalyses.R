#' Summarise biomass, mortality, ANPP and age by target columns
#'
#' @param cohortData `data.table` with columns `c("B", "mortality", "aNPPAct", "age", "noFires", byCols)`
#' @param byCols character. Column names over which summaries are calculated
#'
#' @details Biomass, mortality and ANPP are summed across observations within
#'   combinations of `byCols`, with raw values divided by 100. Age is the mean
#'   biomass-weighted age.
#' @export
makeSummaryTable <- function(cohortData, byCols) {
  summaryCohortData <- cohortData[, list(BiomassBySpecies = as.numeric(sum((B/100), na.rm = TRUE)),
                                         MortalityBySpecies = as.numeric(sum((mortality/100), na.rm = TRUE)),
                                         aNPPBySpecies = as.numeric(sum((aNPPAct/100), na.rm = TRUE)),
                                         avgAgeBySppSim = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                                       sum((B/100), na.rm = TRUE))),
                                  by = byCols]
  summaryCohortData[is.nan(avgAgeBySppSim), avgAgeBySppSim := 0]
  summaryCohortData[, firePresAbs := as.integer(noFires != 0)]
  summaryCohortData
}

#' Prepare data for species-level age comparisons
#'
#' Calculates median ages per species or pixel (stand)
#'
#' @param simData `data.table`. All simulated cohort age data to be considered.
#'   With columns `c("scenario", "rep", "pixelIndex", "year", "firePresAbs", "vegTypeCN", "speciesCode", "B" and "age")`
#' @param obsData `data.table`. All observed tree and stand age data. With columns
#'   `c("ageAtMinDBH", "Reconstructed.age", "Stand.age.3", "speciesCode", "Cover.dendro", "firePresAbs")`
#' @param addLandscape `logical`. If `TRUE`, median ages will be calculated across
#'   the whole landscape.
#' @param speciesLevel `logical`. When `addLandscape == TRUE`, if `speciesLevel`
#'   is `TRUE` a new "landscape" species is added to `speciesCode`, with averages calculated across all species;
#'   otherwise a "landscape" forest type is added to `vegTypeCN` (with averages across forest types).
#' @param ... passed to `reproducible::Cache`
#'
#' @return a `data.table` with columns:
#'   - age: the raw simulated cohort ages per pixel
#'   - ageWeighted: the simulated biomass-weighted species ages per pixel
#'   - sppAgeSim: the simulated *median* species age per pixel
#'   - standAgeSim: the simulated *median* biomass-weighted stand age per pixel
#'   - avgAgeSppObs: the *median* observed species age (basal-area-weighted)
#'     calculated across patches per species
#'   - avgAgeStandObs: the *median* observed stand age calculated across patches
#'     per forest type
#'
#' @details Cohorts with biomass low biomass at old ages and cohorts with ages lower than age at minimum DBH
#'   (`age < ageAtMinDBH`) are excluded from calculations. Low biomass at old age was defined as a biomass
#'   value equal to minimum biomass at `age < ageAtMinDBH`, per species. "Old age" was defined as 0.9*longevity
#'   of the species.
#'   Pixels whose forest type differs from forest types sampled in the field
#'   (`vegTypeCN %in% unique(obsData$Cover.dendro)`) were excluded.
#'   Records in the observed data that correspond to species absent from the simulations
#'   are also excluded.
#'
#'   Species-level median ages (`"sppAgeSim"` and `"avgAgeSppObs"`) are
#'   calculated on raw simulated (`"age"`) and observed (`"Reconstructed.age"` in `obsData`)
#'   ages by species within each simulated pixel or across observed patches, respectively.
#'
#'   Stand-level median ages (`"standAgeSim"` and `"avgAgeStandObs"`) are calculated as
#'   by pixel (`"standAgeSim"`), and across patches (patch.ID) per forest type (`"vegTypeCN"`).
#'   Simulated stand ages (`"standAgeSim"`) are the median biomass-weighted cohort ages in a pixel (with cohort biomass/100 as
#'   weights); median observed stand ages (`"avgAgeStandObs"`) are baed on basal-area-weighted ages in a patch (`"Stand.age.3"` in `obsData`).
#'
#'   All calculations on simulated data are done per year, scenario and replicate.
#'
#' @export
ageComp_data <- function(simData, obsData, speciesTraits, addLandscape = FALSE, speciesLevel = TRUE, ...) {
  ## subset obsData
  obsData <- copy(obsData)
  setnames(obsData, "Cover.dendro", "vegTypeCN", skip_absent = TRUE)
  allCols <- c("Reconstructed.age", "Stand.age.3", "ageAtMinDBH", "speciesCode", "vegTypeCN", "firePresAbs")
  obsData <- obsData[, ..allCols]
  obsData <- obsData[!is.na(speciesCode)]

  ## subset simData
  simData <- simData[B > 0 & vegTypeCN %in% unique(obsData$vegTypeCN)]
  simData <- unique(obsData[, .(speciesCode, ageAtMinDBH)])[simData, on = .(speciesCode)]
  simData <- simData[age >= ageAtMinDBH,]   ## to match field methods

  ## there are some very old cohorts that have virtually no biomass, when compared to very young cohorts
  ## remove cohorts with too low biomass
  minBs <- copy(simData)
  minBs[, minAge := min(age), by = speciesCode]
  minBs <- minBs[age == minAge][, minB := min(B, na.rm = TRUE), by = speciesCode]
  minBs <- unique(minBs[, .(speciesCode, minAge, minB)])
  simData <- minBs[simData, on = "speciesCode"]
  simData <- speciesTraits[, .(speciesCode, longevity)][simData, on = .(speciesCode)]
  simData[, longevity := longevity * 0.9]  ## exclude small cohorts near longevity
  simData[, exclude := FALSE]
  simData[B <= minB & age > longevity, exclude := TRUE] ## there are some very old cohorts that have virtually no biomass
  simData <- simData[exclude == isFALSE(exclude)]
  simData[, `:=`(ageAtMinDBH = NULL,
                 minAge = NULL,
                 minB = NULL,
                 longevity = NULL,
                 exclude = NULL)]

  ## calculate weighted cohort ages first
  ## don't average across years/pixels -- leave all replicates.
  cols <- which(as.matrix(simData[, lapply(.SD, is.numeric)]))
  cols <- names(simData)[cols]
  cacheExtra <- colSums(simData[, ..cols])
  byCols <- c("scenario", "rep", "pixelIndex", "year", "firePresAbs", "vegTypeCN", "speciesCode")
  simData2 <- Cache(.biomassWeightedAge,
                    cohortData = simData,
                    byCols = byCols,
                    .cacheExtra = list(cacheExtra),
                    omitArgs = c("cohortData", "userTags"),
                    userTags = c("biomassWeightedAge", byCols),
                    ...)
  simData <- simData2[simData, on = byCols]

  if (any(is.na(simData$age))) {
    stop("There should be no NAs")
  }

  if (any(is.na(simData$ageWeighted))) {
    stop("There should be no NAs")
  }

  ## repeat values for landscape -- best after slow calculation of biomass-weighted
  ## ages
  if (addLandscape) {
    simDataLandscape <- copy(simData)
    obsDataLandscape <- copy(obsData)

    if (speciesLevel) {
      simDataLandscape[, speciesCode := "landscape"]
      obsDataLandscape[, speciesCode := "landscape"]
    } else {
      simDataLandscape[, vegTypeCN := "landscape"]
      obsDataLandscape[, vegTypeCN := "landscape"]
    }

    simData <- rbind(simData, simDataLandscape)
    obsData <- rbind(obsData, obsDataLandscape)
  }

  ## simulated median biomass-weighted pixel/stand age
  simData <- simData[, standAgeSim := median(ageWeighted, na.rm = TRUE),
                     by = .(scenario, rep, year, pixelIndex)]
  ## simulated median species-age per pixel
  simData <- simData[, sppAgeSim := median(age, na.rm = TRUE),
                     by = .(scenario, rep, year, pixelIndex, speciesCode)]

  ## observed median species-age across patches by species/landscape
  obsData <- obsData[, avgAgeSppObs := median(Reconstructed.age, na.rm = TRUE),
                     by = .(speciesCode, firePresAbs)]
  ## observed median basal-area-weighted (see stand.age.func)
  ## species-age across patches by forest type/landscape
  obsData <- obsData[, avgAgeStandObs := median(Stand.age.3, na.rm = TRUE),
                     by = .(vegTypeCN, firePresAbs)]

  simData <- unique(simData)

  ## join simulated and observed
  onCols <- c("speciesCode", "vegTypeCN", "firePresAbs")
  obsData <- unique(obsData[, .SD, .SDcols = c(onCols, "avgAgeSppObs", "avgAgeStandObs")])

  allData <- obsData[simData, on = onCols]
  allData[, firePresAbs := as.factor(firePresAbs)]
  allData
}


.biomassWeightedAge <- function(cohortData, byCols) {
  cohortData[,  ## absences can't lower the average ages
             list(ageWeighted = as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                             sum((B/100), na.rm = TRUE))),
             by = byCols]
}


#' Summary of vegetation dynamics - plots
summaryPlot <- function(data, x, y, colour, xlabels = NULL, colValues, colLabels,
                        titleLab = "", subtitleLab = NULL, xLab = "", yLab = "",
                        fun.data = "mean_sd", x.text.angle = 0) {
  plotOut <- ggplot(data = data,
                    aes(x = !!sym(x), y = !!sym(y), colour = !!sym(colour))) +
    stat_summary(fun.data = fun.data, position = position_dodge(width = 0.8),
                 linewidth = 1) +
    theme_pubr(base_size = 16, legend = "bottom", x.text.angle = x.text.angle) +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    scale_colour_manual(values = colValues)
  if (!is.null(xlabels)) {
    plotOut <- plotOut +
      if (is(data[[x]], "numeric")) {
        scale_x_manual(labels = xlabels)
      } else {
        scale_x_discrete(labels = xlabels)
      }
  }
  plotOut <- plotOut +
    labs(title = titleLab, y = yLab, x = xLab, subtitle = subtitleLab)
  plotOut
}

#' Simulated and observed distributions of focal variables
simObsDistsPlot <- function(simData, x, ySim, colSim,
                            obsData = NULL, yObs = NULL, xlabels,
                            colValues, colLabels, x.text.angle = 0, showMeans = TRUE,
                            titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  plotOut <- ggplot() +
    geom_violin(data = simData,
                mapping = aes(x = !!sym(x), y = !!sym(ySim),
                              fill = !!sym(colSim), colour = !!sym(colSim)),
                alpha = 0.3, position = position_dodge(width = 0.9))
  if (showMeans) {
    plotOut <- plotOut +
      stat_summary(data = simData,
                   mapping = aes(x = !!sym(x), y = !!sym(ySim), colour = !!sym(colSim)),
                   fun.data = "median_hilow_",
                   position = position_dodge(width = 0.9))
  }

  if (!is.null(obsData)) {
    if (class(simData[[x]]) != class(obsData[[x]])) {
      message("class(simData[[x]]) != class(obsData[[x]]), will try to coherce the later")
      coherceFun <- paste0("as.", class(simData[[x]]), "(obsData[[x]])")
      obsData[[x]] <- eval(parse(text = coherceFun))
    }

    plotOut <- plotOut +
      geom_violin(data = obsData,
                  mapping = aes(x = !!sym(x), y = !!sym(yObs),
                                colour = "observed", fill = "observed"),
                  alpha = 0.3)
    if (showMeans) {
      plotOut <- plotOut +
        stat_summary(data = obsData,
                     aes(x = !!sym(x), y = !!sym(yObs), colour = "observed"),
                     fun.data = "median_hilow_")
    }
  }

  if (!is.null(xlabels)) {
    plotOut <- plotOut +
      if (is(simData[[x]], "numeric")) {
        scale_x_manual(labels = xlabels)
      } else {
        scale_x_discrete(labels = xlabels)
      }
  }

  plotOut <- plotOut +
    scale_colour_manual(values = colValues, labels = colLabels, guide = "none") +
    scale_fill_manual(values = colValues, labels = colLabels) +
    theme_pubr(base_size = 16, legend = "bottom", x.text.angle = x.text.angle) +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab)
  plotOut
}

#' Deviations plots of focal variables from observed values
DevPlot <- function(data, x, y, fill, xlabels, xorder = NULL,
                    fillValues, fillLabels, xFacet, labllr,
                    titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  plotOut <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill))) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey") +
    geom_boxplot() +
    theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +

    if (is.null(xorder)) {
      scale_x_discrete(labels = xlabels)
    } else {
      scale_x_discrete(labels = xlabels, limits = xorder)
    }

  plotOut <- plotOut +
    scale_fill_manual(values = fillValues, labels = fillLabels) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab)
  if (!missing(xFacet)) {
    plotOut <- plotOut +
      facet_grid(rows = xFacet, labeller = labllr)
  }

  plotOut
}

#' MAD plots of focal variables
MADPlot <- function(data, x, y, colour, colValues, colLabels,
                    labllr, xlabels, xorder, xFacet = NULL, x.text.angle = 0,
                    titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  plotOut <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), colour = !!sym(colour), shape = !!sym(colour))) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.5)) +
    scale_y_continuous(limits = c(0, max(data[[y]]))) +
    scale_colour_manual(values = colValues, labels = colLabels, name = "") +
    scale_shape_discrete(labels = colLabels, name = "") +
    scale_x_discrete(labels = xlabels, limits = rev(xorder)) +   ## rev because we flip axes
    theme_pubr(base_size = 16, legend = "bottom", x.text.angle = x.text.angle) +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted"),
          panel.grid.major.x = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab) +
    coord_flip()
  if (!is.null(xFacet)) {
    plotOut <- plotOut +
      facet_grid(rows = xFacet, labeller = labllr)
  }
  plotOut
}


#' Density plots of focal variables
densityPlot <- function(data, x, fill, alpha, fillValues, fillLabels, alphaValues,
                        titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  ggplot(data,
         aes(x = !!sym(x), fill = !!sym(fill), alpha = !!sym(alpha))) +
    geom_density(adjust = 1.5) +
    scale_x_log10() +
    scale_alpha_manual(values = alphaValues, guide = "none") +
    scale_fill_manual(values = fillValues, labels = fillLabels) +
    theme_pubr(base_size = 14, legend = "bottom") +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab)
}


#' Boxplots of focal variables
boxPlot <- function(data, x, y, fill, xlabels, xorder, fillValues, fillLabels,
                    titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  ggplot(data) +
    geom_boxplot(aes(x = "landscape", y = !!sym(y), fill = !!sym(fill)),
                 outlier.size = 0.5) +
    geom_boxplot(aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill)),
                 outlier.size = 0.5) +
    scale_y_log10() +
    scale_x_discrete(labels = xlabels, limits = xorder) +
    scale_fill_manual(values = fillValues, labels = fillLabels) +
    theme_pubr(base_size = 14, legend = "bottom") +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab)
}


#' Species composition plots
sppCompositionPlots <- function(data, x, y, fill, fillValues, fillLabels,
                                xFacet, xlabels, labllr, errorBar = FALSE,
                                titleLab = "", subtitleLab = NULL, xLab = "", yLab = "") {
  plotOut <- ggplot(data, aes(x = !!sym(x), y = !!sym(y), fill = !!sym(fill)))

  if (errorBar) {
    plotOut <- plotOut +
      stat_summary(fun = "mean", geom = "bar", position = position_dodge(width = 0.9)) +
      stat_summary(fun.data = "mean_sd", geom = "linerange", position = position_dodge(width = 0.9))
  } else {
    plotOut <- plotOut +
      stat_summary(fun = "mean", geom = "bar", position = "fill")
  }
  plotOut <- plotOut +
    scale_fill_manual(values = fillValues, labels = fillLabels) +
    scale_x_discrete(labels = xlabels)

  if (!errorBar) {
    plotOut <- plotOut +
      scale_y_continuous(breaks = c(seq(0, 1, 0.25)),
                         labels = c(as.character(seq(0, 1, 0.25))))
  }

  plotOut <- plotOut +
    theme_pubr(base_size = 16, legend = "bottom", x.text.angle = 45) +
    theme(legend.title = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(title = titleLab, y = yLab, subtitle = subtitleLab, x = xLab) +
    facet_grid(rows = xFacet,
               labeller = labllr)

  plotOut
}
#' Hypervolume size boxplots
HVBoxplots <- function(plotData, x = "vegType", y = "logVolume", fill = "scenario", alpha,
                       yLab = y, xLab = x, fillLab = fill,
                       alphaLab, titleLab = "",
                       vegType, HVtype, logVolume, scenario, xLabels, fillLabels, fillVals) {
  plotOut <- ggplot(plotData,
                    aes(x = !!sym(x), y = !!sym(y)
                        # , alpha = scenario
                        , fill = !!sym(fill))) +
    geom_boxplot()

  if (!missing(xLabels)) {
    plotOut <- plotOut +
      scale_x_discrete(labels = xLabels)
  }

  if ((!missing(fillLabels) & missing(fillVals)) |
      missing(fillLabels) &!missing(fillVals)) {
    stop("fillLabels and fillVals must be both provided, or both not provided.")
  }

  if (!missing(fillLabels) & !missing(fillVals)) {
    plotOut <- plotOut +
      scale_fill_manual(labels = fillLabels, values = fillVals)
  }

  if (!missing(alpha)) {
    plotOut <- plotOut +
      scale_alpha_manual(values = c("noPM" = 0.4, "PM" = 1.0), labels = scenLabels)
    if (missing(alphaLab))
      alphaLab <- alpha
  }

  if (missing(alphaLab))
    alphaLab <- "alpha"

  plotOut <- plotOut +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    theme_pubr(base_size = 12, margin = FALSE) +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(x = xLab, y = yLab, fill = fillLab, alpha = alphaLab, title = titleLab)
  # guides(alpha = guide_legend(override.aes = list(fill = "grey50"))) +
  # facet_wrap(~ HVtype, nrow = 2, scales = "free_y",
  #            labeller = labeller(HVtype = c("vegHV" = "Forest diversity",
  #                                           "fireHV" = "Pyrodiversity")))

  plotOut
}
#' Biodiversity ~ Pyrodiversity relationships plots
#' I.e. veg HV size ~ fire HV size
#' ... is passed to stat_smooth
plotBioPyroFunSmooth <- function(plotData, x = "logFireHV", yPoints = "logVegHV",
                                 linetype = "scenario", colour = "vegType",
                                 yPred = "pred", ymin, ymax,
                                 #shape = "scenario"
                                 colourLabels, colourVals,
                                 titleLab = "", xLab = "Pyrodiversity", yLab = "Forest diversity",
                                 colourLab = "", ...) {
  # old code, but here to remind of different formulae
  # if (all(plotData$vegType == "landscape")) {
  #   form <- quote(y ~ x)
  # } else {
  # form <- quote(y ~ x + I(x^2))   ## needs x = logFireHVcenter below
  # form <- quote(y ~ poly(x, 2))
  # }

  plotOut <- ggplot(plotData,
         aes(x = !!sym(x), linetype = !!sym(linetype), colour = !!sym(colour)
           #, shape = !!sym(shape)
         )) +
    geom_point(aes(y = !!sym(yPoints)))

  xLimits <- range(plotData[[x]])

  plotOut <- plotOut +
    stat_smooth(aes(y = !!sym(yPred)), ...)

  if (!missing(ymin)) {
    if (missing(ymax)) {
      stop("Provide both ymin and ymax")
    }
    plotOut <- plotOut +
      geom_ribbon(aes(ymin = !!sym(ymin), ymax = !!sym(ymax), fill = !!sym(colour)),
                  alpha = 0.5)
  }

  if (!missing(colourLabels)) {
    if (missing(colourVals)) {
      stop("Provide both colourLabels and colourVals")
    }
    plotOut <- plotOut +
      scale_colour_manual(labels = colourLabels, values = colourVals)
  }

  plotOut <- plotOut +
    # scale_linetype_manual(labels = scenLabels,
    #                       values = scenLinetype) +
    # scale_shape_discrete(labels = scenLabels) +
    scale_x_continuous(limits = xLimits) +
    theme_pubr(base_size = 12, margin = TRUE) +
    theme(legend.box = "vertical",
          strip.background = element_blank(),
          panel.grid.major.y = element_line(colour = "grey", linewidth = 11/22, linetype = "dotted")) +
    labs(x = xLab, y = yLab, title = titleLab, colour = colourLab
         #linetype = linetypeLab, shape = shapeLab
    ) +
    # facet_wrap( ~ vegType, labeller = labeller(vegType = vegTypeCNLabels),
    #             scales = "free") +
    guides(linetype = guide_legend(override.aes = list(colour = "black")))

  plotOut
}

#' Summarise fire regime attributes per pixel
#'
#' As in Steel et al 2021, fire properties are summarized across time, but by pixel
#' (and by scenario/rep) as pixel-level averages.
#'
#' Only pixels with vegetation dynamics are considered so that we can compare with biodiv. HVs
#'
#' @param allPixelBurnData `data.table` with all year and pixel-level fire regime attributes:
#'    -  `fireFreq` (fire frequency) -- already temporally integrated so the "summary" is the
#'       unique value per pixel.
#'    -  `severity` (fire severity in category)
#'    -  `severityB` (fire severity as B lost)
#'    -  `patchSizeLogHa` (severity patch size in log-Ha units)
#' @param allPixelCohortDataMnt `data.table` with all year and pixel-level cohort B and age.
#'    Only used to filter `allPixelBurnData`
#' @param addMedian logical. Should medians be computed?
#'    -  `meanFreq` / `medianFreq` (the unique value of `fireFreq`)
#'       unique value per pixel.
#'    -  `meanSev` / `medianSev` (mean/median `severity`)
#'    -  `meanSevB` / `medianSevB` (mean/median `severityB`)
#'    -  `meanPatchS` / `medianPatchS` (mean/median `patchSizeLogHa`)
#' @returns a table with the four fire regime properties averaged across time and by pixel:
#'
#' @import data.table
summariseFireRegAttrs <- function(allPixelBurnData, allPixelCohortDataMnt, yearSubset, addMedian = FALSE) {

  summaryFireAttributes <- allPixelBurnData[pixelIndex %in% allPixelCohortDataMnt$pixelIndex]

  if (addMedian) {
    summaryFireAttributes <- summaryFireAttributes[, list(meanFreq = unique(fireFreq),   ## note that fireFreq is already an average of mean fire intervals
                                                          meanSev = mean(severity),
                                                          meanSevB = mean(severityB),
                                                          meanSevPropB = mean(severityPropB),
                                                          meanPatchS = mean(patchSizeLogHa),
                                                          medianFreq = unique(fireFreq),   ## note that fireFreq is already an average of median fire intervals
                                                          medianSev = median(severity),
                                                          medianSevB = median(severityB),
                                                          medianSevPropB = median(severityPropB),
                                                          medianPatchS = median(patchSizeLogHa)),
                                                   by = .(scenario, rep, pixelIndex)]
  } else {
    summaryFireAttributes <- summaryFireAttributes[, list(meanFreq = unique(fireFreq),   ## note that fireFreq is already an average of mean fire intervals
                                                          meanSev = mean(severity),
                                                          meanSevB = mean(severityB),
                                                          meanSevPropB = mean(severityPropB),
                                                          meanPatchS = mean(patchSizeLogHa)),
                                                   by = .(scenario, rep, pixelIndex)]
  }


  ## add vegType per pixel at the end of the simulation
  ## and add pixels that had no fires
  cols <- c("pixelIndex", "vegTypeCN", "scenario", "rep")
  summaryFireAttributes <- summaryFireAttributes[unique(allPixelCohortDataMnt[year == max(yearSubset), ..cols]),
                                                 on = c("scenario", "rep", "pixelIndex")]
  ## checks
  if (getOption("LandR.assertions", TRUE)) {
    test1 <- any(is.na(summaryFireAttributes$vegTypeCN))
    if (test1) {
      stop("NA vegTypeCNs where pixelGroup - i.e. vegetation - exists")
    }

    test2 <- allPixelBurnData[pixelIndex %in% allPixelCohortDataMnt$pixelIndex][summaryFireAttributes[is.na(meanFreq), .(scenario, rep, pixelIndex)],
                                                                                on = .(scenario, rep, pixelIndex), nomatch = 0]
    if (nrow(test2)) {
      stop("pixels that had fire and veg data in allPixelBurnData were ",
           "accidentally dropped when adding vegTypeCN")
    }

    test3 <- any(is.na(summaryFireAttributes$meanFreq))
    test4 <- any(is.na(summaryFireAttributes$meanSev))
    test5 <- any(is.na(summaryFireAttributes$meanSevB))
    test6 <- any(is.na(summaryFireAttributes$meanPatchS))

    if (any(test3, test4, test5, test6)) {
      stop("Found NAs in fire properties")
    }

    rm(test1, test2, test3, test4, test5, test6)
  }
  gc(reset = TRUE)

  return(summaryFireAttributes)
}


#' VEGETATION DATA FOR HVs
#'
#' Prepares vegetation data (relative tree species biomass, stand biomass,
#' mean and SD of stand age) for hypervolume calculations.
#'
#' @param allPixelCohortDataMnt
#' @param summaryFireAttributes
#' @param useFirstLastYear logical. Should only the first and last
#'   years of the simulation be used?
#' @param yearSubset Numeric vector. If `useFirstLastYear == TRUE`, this is the
#'   a vector of simulation years used to run analyses, from which the minimum (first) and
#'   maximum (last) years will be subset to prepare the vegetation data.
#' @param yearSamples Optional. Numeric vector. If `useFirstLastYear == FALSE` and
#'   Sample of simulation years to be used to prepare the vegetation data. If not
#'   If not provided, all years from `yearSubset` will be used.
#'
#' @details Mean and SD of stand age (mean/sdStandAge) are the mean and SD of biomass-weighted
#'   cohort ages in a pixel (calculated across species). NA values are converted to
#'   0s, as they indicate lack of forest cover.
#'   NaN relative biomasses (relB) are converted to 0, as they indicate lack of forest cover in
#'   the pixel.
#'   All output vegetation attributes (relative species biomass,mean/SD of stand age)
#'   are calculated per pixel, year, repetition and scenario
#'
#' @returns a `data.table` with columns:
#'   * scenario, rep, year, pixelIndex, speciesCode, vegTypeCN -- identifiers
#'   * meanStandAge, sdStandAge, relB
#' @export
#' @import data.table
prepVegDataHVs <- function(allPixelCohortDataMnt, summaryFireAttributes,
                           useFirstLastYear, yearSubset, yearSamples) {

  if (useFirstLastYear) {
    vegDataForHVs <- allPixelCohortDataMnt[year %in% c(min(yearSubset), max(yearSubset))]
  } else {
    if (exists("yearSamples")) {
      vegDataForHVs <- allPixelCohortDataMnt[yearSamples, on = .(year, rep)]
    } else {
      vegDataForHVs <- allPixelCohortDataMnt[year %in% yearSubset]
    }
  }

  if (getOption("LandR.assertions", TRUE)) {
    pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
    temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
    setkey(temp, scenario, rep, pixelIndex)
    setkey(vegDataForHVs, scenario, rep, pixelIndex)

    if (isFALSE(identical(temp, vegDataForHVs))) {
      stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
           "Combinations as allPixelCohortDataMnt")
    }

    temp <- vegDataForHVs[, length(unique(pixelIndex)), by = .(scenario, rep, year)]
    if (length(unique(temp$V1)) > 1) {
      stop("There should be the same number of pixels every year.")
    }

    temp <- split(vegDataForHVs[, .(pixelIndex, scenario, rep, year)],
                  by = c("scenario", "rep", "year"), keep.by = FALSE)
    temp <- lapply(temp, FUN = function(x) unique(x[["pixelIndex"]]))
    test <- lapply(1:length(temp), function(n) setdiff(temp[[n]], unlist(temp[-n])))

    test <- sapply(test, length)
    if (any(test))
      stop("Different pixelIndex between scenario/rep/year combinations")

    if (min(vegDataForHVs$year) == 2011) {
      temp <- split(vegDataForHVs[year == min(yearSubset), .(vegTypeCN, pixelIndex, scenario, rep)],
                    by = c("scenario", "rep"), keep.by = FALSE)
      temp <- lapply(temp, FUN = function(x) setkey(x, vegTypeCN, pixelIndex))
      temppix <- lapply(temp, FUN = function(x) x[["pixelIndex"]])
      tempveg <- lapply(temp, FUN = function(x) x[["vegTypeCN"]])

      test <- lapply(1:length(temppix), function(n) setdiff(temppix[[n]], unlist(temppix[-n])))
      test2 <- lapply(1:length(tempveg), function(n) setdiff(tempveg[[n]], unlist(tempveg[-n])))
      test <- sapply(test, length)
      test2 <- sapply(test2, length)

      if (any(test) | any(test2))
        stop("Difference pixelIndex/vegTypeCN combinations between scenario/reps in the first year")
    }

    ## checks at landscape scale:
    pixelIndices <- unique(summaryFireAttributes[,.(scenario, rep, pixelIndex)])
    temp <- vegDataForHVs[pixelIndices, on = .(scenario, rep, pixelIndex), nomatch = 0]
    setkey(temp, scenario, rep, pixelIndex)
    setkey(vegDataForHVs, scenario, rep, pixelIndex)

    if (isTRUE(any(temp != vegDataForHVs))) {
      stop("Something is wrong. summaryFireAttributes should have the same pixelIndex/scenario/rep\n",
           "Combinations as allPixelCohortDataMnt")
    }

    rm(temp, test, test2)
    for(i in 1:3) gc(reset = TRUE)
  }

  ## prep data for hypervolumes
  ## calculate stand-age as the mean biomass-weighted age
  vegDataForHVs[, `:=`(meanStandAge = mean(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                        sum((B/100), na.rm = TRUE))),
                       sdStandAge = sd(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                    sum((B/100), na.rm = TRUE)))),
                by = .(scenario, rep, year, pixelIndex, vegTypeCN)]
  vegDataForHVs[is.na(meanStandAge), meanStandAge := 0]  ## NAs come from stands with 0 B and 0 age
  vegDataForHVs[is.na(sdStandAge), sdStandAge := 0]  ## NAs come from stands with 0 B and 0 age or just one value of standAge

  ## calculate relative species B (across speciesCohorts)
  vegDataForHVs[, standB := sum(B, na.rm = TRUE), by =  c("scenario", "rep", "year", "pixelIndex", "vegTypeCN")]
  vegDataForHVs[, relB := sum(B) / standB, by = c("scenario", "rep", "year", "pixelIndex", "vegTypeCN", "speciesCode")]
  vegDataForHVs[relB == "NaN", relB := 0]

  ## expand data
  cols <- c("meanStandAge", "sdStandAge", "relB", "speciesCode", "scenario", "rep", "year", "pixelIndex", "vegTypeCN")   ## keep rep for wrapper.
  vegDataForHVs <- unique(vegDataForHVs[, ..cols])
  vegDataForHVs <- dcast.data.table(vegDataForHVs, as.formula("... ~ speciesCode"),
                                    value.var = "relB")

  return(vegDataForHVs)
}
