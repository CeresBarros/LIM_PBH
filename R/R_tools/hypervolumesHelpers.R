#' Wrapper function to calculate vegetation attributes hypervolumes
vegHVWrapper <- function(allData, HVcols, IDcols, HVIDcol, file.suffix, addNoise = TRUE, ...) {
  ## a bit of prep
  HVnames <- unique(allData[[HVIDcol]])
  if (length(HVnames) != 2) {
    stop("There should be 2 values in,", HVIDcol)
  }

  ## add noise to data if all dimensions have zero variance
  if (addNoise) {
    IDcols2 <- unique(c(HVIDcol, setdiff(IDcols, "pixelIndex")))
    needsNoise <- allData[, lapply(.SD, sd), .SDcols = HVcols, by = IDcols2]
    needsNoise <- needsNoise[,  rowSums(.SD, na.rm = TRUE) == 0, .SDcols = HVcols, by = IDcols2]

    if (any(needsNoise$V1)) {
      tempData <- allData[needsNoise[V1 == TRUE, ..IDcols2],
                          on = IDcols2]
      tempData[, (HVcols) := lapply(.SD, function(x, n) rnorm(n, mean(x), 0.0000001),
                                    n = .N), .SDcols = HVcols, by = IDcols2]
      allData <- rbind(allData[!tempData, on = IDcols2],
                       tempData, fill = TRUE, use.names = TRUE)
    }
  }

  ## HV calculation
  out <- tryCatch(hypervolumes(HVdata1 = as.data.frame(allData[get(HVIDcol) == HVnames[1]]),
                               HVdata2 = as.data.frame(allData[get(HVIDcol) == HVnames[2]]),
                               HVidvar = which(names(allData) == HVIDcol),
                               init.vars = which(names(allData) %in% HVcols),
                               file.suffix = file.suffix,
                               ...), error = function(e) e)

  if (is(out, "error")) {
    if (grepl("contour.default", out)) {
      hypervolumes(HVdata1 = as.data.frame(allData[get(HVIDcol) == HVnames[1]]),
                   HVdata2 = as.data.frame(allData[get(HVIDcol) == HVnames[2]]),
                   HVidvar = which(names(allData) == HVIDcol),
                   init.vars = which(names(allData) %in% HVcols),
                   file.suffix = file.suffix,
                   plotHVDots = list(contour.type = "ball"),
                   ...)
    }
  }
}

#' Wrapper function to calculate fire attributes hypervolumes
fireHVWrapper <- function(allData, cols, file.suffix, addNoise = TRUE, ...) {
  print(file.suffix)

  if (addNoise) {
    ## add noise to data if all dimensions have zero variance
    needsNoise <- allData[, lapply(.SD, sd), .SDcols = cols, by = .(scenario, rep)]
    needsNoise <- needsNoise[,  rowSums(.SD, na.rm = TRUE) == 0, .SDcols = cols, by = .(scenario, rep)]

    if (any(needsNoise$V1)) {
      tempData <- allData[needsNoise[V1 == TRUE, .(scenario, rep)],
                          on = .(scenario, rep)]
      tempData[, (cols) := lapply(.SD, function(x, n) rnorm(n, mean(x), 0.0000001),
                                  n = .N), .SDcols = cols]
      allData <- rbind(allData[!tempData, on = .(scenario, rep)],
                       tempData, fill = TRUE, use.names = TRUE)
    }
  }

  out <- tryCatch(hypervolumes(HVdata1 = as.data.frame(allData[scenario == "noPM"]),
                               HVdata2 = as.data.frame(allData[scenario == "PM"]),
                               HVidvar = which(names(allData) == "scenario"),
                               init.vars = which(names(allData) %in% cols),
                               file.suffix = file.suffix,
                               ...), error = function(e) e)

  if (is(out, "error")) {
    if (grepl("contour.default", out)) {
      hypervolumes(HVdata1 = as.data.frame(allData[scenario == "noPM"]),
                   HVdata2 = as.data.frame(allData[scenario == "PM"]),
                   HVidvar = which(names(allData) == "scenario"),
                   init.vars = which(names(allData) %in% cols),
                   file.suffix = file.suffix,
                   plotHVDots = list(contour.type = "ball"),
                   ...)
    }
  }
}


#' Function to load hypervolume comparison tables from .rds files that match a pattern
#'
#' @param x the pattern of file name to be matched
#' @param files the vector of all file names to be searched

loadHVResultsFromRDS <- function(x, files) {
  files <- grep(x, files, value = TRUE)

  HVData <- lapply(files, FUN = function(ff) {
    HVData <- readRDS(ff)
    if (!is(HVData, "data.table")) {
      ## intersection results only have one row
      if (!is.null(rownames(HVData)) & nrow(HVData) > 1) {
        HVData$HVid <- rownames(HVData)
      }
      HVData <- as.data.table(HVData)
    }
    ff <- basename(ff)
    if (grepl("yr", ff)) {
      yr <- sub(".*yr([0-9]{1,4})_.*", "\\1", ff)
    } else {
      yr <- NA
    }
    if (grepl("rep", ff)) {
      r <- sub(".*(rep)([0-9]+)_.*", "\\2", ff)
    } else {
      r <- NA
    }

    if (grepl("[0-9]+\\.rds$", ff)) {
      rHV <- sub(".*([0-9]+)\\.rds", "\\1", ff)
    } else {
      r <- NA
    }

    vegType <- unlist(strsplit(sub("(_yr|_rep).*", "", ff), "_"))[2]
    scenario <- unlist(strsplit(sub("(_yr|_rep).*", "", ff), "_"))[3]

    HVData[, year := as.integer(yr)]
    HVData[, rep := as.integer(r)]
    HVData[, repHV := as.integer(rHV)]
    HVData[, vegType := vegType]
    HVData[, scenario := scenario]

    return(HVData)
  }) %>%
    rbindlist(fill = TRUE, l = ., use.names = TRUE)
  return(HVData)
}


#' Wrapper function to plot hypervolumes with PCA factor vectors
#' and vectors fitted post-hoc
#'
#' @param vegType
#' @param vegHVPCAscores a table of factor scores from PCA
#' @param pixelIndexDT a data.table with pixelIndices to track across hypervolumes
#'   (should be unique per combination of repetition and vegetation type)
#' @param vegTypeCNLabels named vector containing vegType labels (names must correspond to \code{vegType}.
#' @param mergeVegType character. Can be "mergeDMCPSME", "mergePSME" of NULL to merge
#'   dry mixed conifer and any Doug-fir dominated stands, just Doug-fir dominated stands
#'   or no merging.
#' @param colsHV names of PCA axes used to make hypervolumes
#'   Used to for plot title
#' @param cacheRepo passed to \code{reproducible::Cache}
#' @param figOutputPath folder where to save hypervolume plot.
#' @param ... more arguments passed to \code{ToolsCB::plotHypervolumes3D}

plotHVs3DWrapper <- function(vegType, vegHVPCAscores, pixelIndexDT, vegTypeCNLabels,
                             mergeVegType = NULL, startYear, endYear, colsHV, cacheRepo, figOutputPath,
                             ...) {
  if (is(vegType, "factor")) {
    vegType <- as.character(vegType)
  }
  grepStr <- vegType

  if (grepl("PSME", vegType)) {
    if (mergeVegType == "mergeDMCPSME") {
      grepStr <- "PSME"
    } else {
      if (mergeVegType == "mergePSME" & grepl("^dryPSME|^PSME", vegType)) {
        grepStr <- "^dryPSME|^PSME"
      }
    }
  }

  HVcombos <- if (is.null(startYear) & is.null(endYear)) {
    expand.grid(c("noPM", "PM"), c(NA))
  } else {
    expand.grid(c("noPM", "PM"), c(startYear, endYear))
  }
  names(HVcombos) <- c("scenario", "year")

  allHVs <- Map(scen = HVcombos$scenario,
                yearToMatch = HVcombos$year,
                f = .HVforPlots,
                MoreArgs = list(pixelIndexDT = pixelIndexDT,
                                vegHVPCAscores = vegHVPCAscores,
                                colsHV = colsHV,
                                grepStrVegType = grepStr,
                                userTags = vegType,
                                cacheRepo = cacheRepo))

  ## rename using scenario only, for plotting
  allHVs <- lapply(allHVs, function(HV) {
    HV@Name <- sub("_", "", sub(vegType, "", HV@Name))
    HV
    })

  HVls <- hypervolume::hypervolume_join(allHVs)
  args <- c(HVlist = HVls, list(...))
  args$main <- vegTypeCNLabels[vegType]
  names(args$colors) <- sapply(allHVs, function(HV) {HV@Name})
  names(args$centroid.cols) <- names(args$colors)

  if (is.null(args$limits)) {
    HVpoints <- rbindlist(lapply(allHVs, function(HV) as.data.table(HV@RandomPoints)))
    HVpoints <- HVpoints[, 1:3] ## only three first axes are plotted
    allPoints <- rbind(HVpoints,
                      args$loadings_coords[, names(HVpoints)],
                      args$PHvect_coords[, names(HVpoints)])
    limits <- cbind(as.matrix(apply(allPoints, 2, min)),
                    as.matrix(apply(allPoints, 2, max)))
    colnames(limits) <- c("min", "max")

    ## convert to list of two element vectors and round
    args$limits <- apply(limits, 1, function(X) X, simplify = "list")
  }

  png(file.path(figOutputPath, paste0("HV3DplotWvectors_", vegType, ".png")), width = 7, heigh = 7,
       units = "in", res = 300)
  do.call(plotHypervolumes3D, args)
  dev.off()
}


.HVforPlots <- function(pixelIndexDT, vegHVPCAscores, scen, colsHV,
                        grepStrVegType, yearToMatch, userTags, cacheRepo) {
  ## subset tables to correct year, veg type and scenario
  tempPixID <- pixelIndexDT[grepl(grepStrVegType, vegTypeCN)]
  tempPixID <- tempPixID[scenario == scen]
  if (is.na(yearToMatch)) {
    tempData <- vegHVPCAscores[scenario == scen,]
    yearToMatch <- NULL
  } else {
    yearToMatch <- ifelse(is.character(vegHVPCAscores$year), as.character(yearToMatch), as.integer(yearToMatch))
    tempData <- vegHVPCAscores[year %in% as.character(yearToMatch) & scenario == scen,]
  }
  ## subset pixels
  tempData <- tempData[tempPixID, on = .(rep, vegTypeCN, pixelIndex)]

  if (!is.null(userTags)) {
    userTags <- c(userTags, "hypervolume", yearToMatch, scen)
  } else {
    userTags <- c("hypervolume", paste(unique(tempData$vegTypeCN), collapse = "_"), yearToMatch, scen)
  }

  cacheObj <- digest::digest(tempData, algo = "xxhash64")

  HV <- Cache(hypervolume::hypervolume,
              data = tempData[, ..colsHV],
              name = paste(scen, unique(as.character(tempData$vegTypeCN)), sep = "_"),
              method = "svm",
              svm.gamma = 0.01,
              .cacheExtra = cacheObj,
              cacheRepo = cacheRepo,
              userTags = userTags,
              omitArgs = c("data"))
  return(HV)
}



## HYPERVOLUMES BY VEGETATION TYPE -----------------------------------
## BANDWITH ESTIMATES ----
# amc::.gc()
# parallel_wrapper <- function(ncores, summaryFireAttributes, byVars, bw.outputPath) {
#   if (!dir.exists(bw.outputPath)) dir.create(bw.outputPath)
#
#   if (.Platform$OS.type == "windows") {
#     plan("multisession", workers = ncores)
#   } else {
#     plan("multicore", workers = ncores)
#   }
#   bw_estimates <- future_lapply(X = split(summaryFireAttributes, by = byVars),
#                                 FUN = function(DT, bw.outputPath) {
#                                   r <- unique(DT$rep)
#                                   veg <- unique(DT$vegTypeCN)
#                                   message(paste("Calculating PCAs and estimating BWs for: rep", r, "and", veg))
#                                   file.suffix <- paste0("fireHVs_freeBW_", veg, "_rep", r)
#
#                                   init.vars <- grep("mean", names(DT))
#
#                                   DT <- ToolsCB:::.scaleVars(DT, init.vars)
#
#                                   out <- estimateBW_wrapper(as.data.frame(DT),
#                                                             init.vars = init.vars,
#                                                             HVidvar = which(names(DT) == "scenario"),
#                                                             noAxes = 4,
#                                                             ordination = "PCA",
#                                                             file.suffix = file.suffix,
#                                                             outputs.dir = bw.outputPath)
#                                   return(out)
#                                 },
#                                 bw.outputPath = bw.outputPath)
#   future:::ClusterRegistry("stop")
#   return(bw_estimates)
# }
#
# bw_estimates <- Cache(parallel_wrapper,
#                       summaryFireAttributes = summaryFireAttributes,
#                       ncores = 10,
#                       byVars = c("rep", "vegTypeCN"),
#                       bw.outputPath = bw.outputPath,
#                       cacheRepo = cacheRepo,
#                       userTags = c("bw_estimates", "hypervolumes", "vegTypeCN"),
#                       omitArgs = c("userTags", "ncores", "bw.outputPath"))
#
# bw_estimates <- do.call(rbind.data.frame, bw_estimates)
# bw_estimates$vegType <- sub("[[:digit:]]*\\.", "", sub("\\.PC.*", "", row.names(bw_estimates)))
# bw_estimates$rep <- sub("\\..*", "", sub("\\.PC.*", "", row.names(bw_estimates)))
# bw_estimates <- as.data.table(bw_estimates)
# saveRDS(bw_estimates, file.path(bw.outputPath, "BW_estimates_vegType.rds"))
#
# summaryBW <- bw_estimates[, list(SilvermanMean = mean(c(SilvBW_HV1, SilvBW_HV2)),
#                                  StDevMean = mean(c(stdev_HV1, stdev_HV2)),
#                                  SilvermanMax = max(c(SilvBW_HV1, SilvBW_HV2)),
#                                  StDevMax = max(c(stdev_HV1, stdev_HV2))),
#                           by = "PC"]
#
# saveRDS(summaryBW, file.path(bw.outputPath, "BW_MeanMax_vegType.rds"))


## fix bandwidth to max of estimated BW
# summaryBW <- readRDS(file.path(bw.outputPath, "BW_MeanMax_vegType.rds"))
# bwHV <- summaryBW[["StDevMax"]]


## HYPERVOLUMES ACROSS THE LANDSCAPE - only montane belt ----------------
## BANDWITH ESTIMATES ----
# amc::.gc()
# parallel_wrapper <- function(ncores, summaryFireAttributes, byVars, bw.outputPath) {
#   if (!dir.exists(bw.outputPath)) dir.create(bw.outputPath)
#
#   if (.Platform$OS.type == "windows") {
#     plan("multisession", workers = ncores)
#   } else {
#     plan("multicore", workers = ncores)
#   }
#   bw_estimates <- future_lapply(X = split(summaryFireAttributes, by = byVars),
#                                 FUN = function(DT, bw.outputPath) {
#                                   r <- unique(DT$rep)
#                                   message(paste("Calculating PCAs and estimating BWs for: rep", r))
#                                   file.suffix <- paste0("fireHVs_freeBW_landscape_rep", r)
#
#                                   init.vars <- grep("mean", names(DT))
#
#                                   DT <- ToolsCB:::.scaleVars(DT, init.vars)
#
#                                   out <- estimateBW_wrapper(as.data.frame(DT),
#                                                             init.vars = init.vars,
#                                                             HVidvar = which(names(DT) == "scenario"),
#                                                             noAxes = 4,
#                                                             ordination = "PCA",
#                                                             file.suffix = file.suffix,
#                                                             outputs.dir = bw.outputPath)
#                                   return(out)
#                                 },
#                                 bw.outputPath = bw.outputPath)
#   future:::ClusterRegistry("stop")
#   return(bw_estimates)
# }
# bw_estimates <- Cache(parallel_wrapper,
#                       ncores = 10,
#                       summaryFireAttributes = summaryFireAttributes,
#                       byVars = c("rep"),
#                       bw.outputPath = bw.outputPath,
#                       cacheRepo = cacheRepo,
#                       userTags = c("bw_estimates", "hypervolumes", "landscape"),
#                       omitArgs = c("userTags", "ncores", "bw.outputPath"))
#
# bw_estimates <- do.call(rbind.data.frame, bw_estimates)
# bw_estimates$rep <- sub("\\..*", "", sub("\\.PC.*", "", row.names(bw_estimates)))
# bw_estimates <- as.data.table(bw_estimates)
# saveRDS(bw_estimates, file.path(bw.outputPath, "BW_estimates_landscape.rds"))
#
# summaryBW <- bw_estimates[, list(SilvermanMean = mean(c(SilvBW_HV1, SilvBW_HV2)),
#                                  StDevMean = mean(c(stdev_HV1, stdev_HV2)),
#                                  SilvermanMax = max(c(SilvBW_HV1, SilvBW_HV2)),
#                                  StDevMax = max(c(stdev_HV1, stdev_HV2))),
#                           by = "PC"]
#
# saveRDS(summaryBW, file.path(bw.outputPath, "BW_MeanMax_landscape.rds"))
#
# ## fix bandwidth to max of estimated BW
# summaryBW <- readRDS(file.path(bw.outputPath, "BW_MeanMax_landscape.rds"))
# bwHV <- summaryBW[["StDevMax"]]
