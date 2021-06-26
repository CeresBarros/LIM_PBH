#' Wrapper function to calculate hypervolumes
hypervolumesWrapper <- function(allData, noAxes, cols, bwVal1, bwVal2, file.suffix) {
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
  hypervolumes(HVdata1 = as.data.frame(allData[scenario == "noPM"]),
               HVdata2 = as.data.frame(allData[scenario == "PM"]),
               HVidvar = which(names(allData) == "scenario"),
               init.vars = which(names(allData) %in% cols),
               HVmethod = "svm", no.runs = 3,
               # freeBW = FALSE, bwHV1 = bwHV, bwHV2 = bwHV,
               svm.gamma = 0.01,
               do.scale = TRUE,
               noAxes = noAxes, outputs.dir = HVoutputPath,
               file.suffix = file.suffix,
               saveOrdi = TRUE, plotOrdi = TRUE, plotHV = TRUE)
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
