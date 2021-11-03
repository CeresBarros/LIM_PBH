Require("dplyr")

#' Wrapper function to calculate vegetation attributess hypervolumes
vegHVWrapper <- function(allData, IDcols, HVIDcol, file.suffix, ...) {
  ## a bit of prep
  HVnames <- unique(allData[[HVIDcol]])
  if (length(HVnames) != 2) {
    stop("There should be 2 values in,", HVIDcol)
  }

  ## calculate stand-age as the mean biomass-weighted age
  allData[, `:=`(meanStandAge = mean(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                          sum((B/100), na.rm = TRUE))),
                 sdStandAge = sd(as.numeric(sum(age * (B/100), na.rm = TRUE) /
                                                  sum((B/100), na.rm = TRUE)))),
          by = IDcols]
  allData[is.na(meanStandAge), meanStandAge := 0]  ## NAs come from stands with 0 B and 0 age
  allData[is.na(sdStandAge), sdStandAge := 0]  ## NAs come from stands with 0 B and 0 age or just one value of standAge

  ## calculate relative species B (across speciesCohorts)
  allData[, standB := sum(B, na.rm = TRUE), by = IDcols]
  allData[, relB := sum(B) / standB, by = c(IDcols, "speciesCode")]
  allData[relB == "NaN", relB := 0]

  ## expand data
  cols <- c("meanStandAge", "sdStandAge", "relB", "speciesCode", IDcols)   ## keep rep for wrapper.
  allData <- unique(allData[, ..cols])
  allData <- dcast.data.table(allData, as.formula("... ~ speciesCode"),
                              value.var = "relB")

  ## still use 4 axes.
  cols <- setdiff(names(allData), IDcols)

  ## add noise to data if all dimensions have zero variance
  IDcols2 <- setdiff(IDcols, "pixelIndex")
  needsNoise <- allData[, lapply(.SD, sd), .SDcols = cols, by = IDcols2]
  needsNoise <- needsNoise[,  rowSums(.SD) == 0, .SDcols = cols, by = IDcols2]

  if (any(needsNoise$V1)) {
    tempData <- allData[needsNoise[V1 == TRUE, ..IDcols2],
                        on = IDcols2]
    tempData[, (cols) := lapply(.SD, function(x, n) rnorm(n, mean(x), 0.0000001),
                                n = .N), .SDcols = cols, by = IDcols2]
    allData <- rbind(allData[!tempData, on = IDcols2],
                     tempData, fill = TRUE, use.names = TRUE)
  }
  ## HV calculation
  out <- tryCatch(hypervolumes(HVdata1 = as.data.frame(allData[get(HVIDcol) == HVnames[1]]),
                               HVdata2 = as.data.frame(allData[get(HVIDcol) == HVnames[2]]),
                               HVidvar = which(names(allData) == HVIDcol),
                               init.vars = which(names(allData) %in% cols),
                               file.suffix = file.suffix,
                               ...), error = function(e) e)

  if (is(out, "error")) {
    if (grepl("contour.default", out)) {
      hypervolumes(HVdata1 = as.data.frame(allData[get(HVIDcol) == HVnames[1]]),
                   HVdata2 = as.data.frame(allData[get(HVIDcol) == HVnames[2]]),
                   HVidvar = which(names(allData) == HVIDcol),
                   init.vars = which(names(allData) %in% cols),
                   file.suffix = file.suffix,
                   plotHVDots = list(contour.type = "ball"),
                   ...)
    }
  }
}

#' Wrapper function to calculate fire attributes hypervolumes
fireHVWrapper <- function(allData, cols, file.suffix, ...) {
  print(file.suffix)
  ## add noise to data if all dimensions have zero variance
  needsNoise <- allData[, lapply(.SD, sd), .SDcols = cols, by = .(scenario, rep)]
  needsNoise <- needsNoise[,  rowSums(.SD) == 0, .SDcols = cols, by = .(scenario, rep)]

  if (any(needsNoise$V1)) {
    tempData <- allData[needsNoise[V1 == TRUE, .(scenario, rep)],
                        on = .(scenario, rep)]
    tempData[, (cols) := lapply(.SD, function(x, n) rnorm(n, mean(x), 0.0000001),
                                n = .N), .SDcols = cols]
    allData <- rbind(allData[!tempData, on = .(scenario, rep)],
                     tempData, fill = TRUE, use.names = TRUE)
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
#                       cacheRepo = simPaths$cachePath,
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
#                       cacheRepo = simPaths$cachePath,
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
