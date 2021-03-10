## DIAGNOSE MODELBIOMASS ---------------------------------
## allFit is taking too long (>12h and didn't even with first optimizer)
if (FALSE) {
  library(lme4)
  modBiomass <- LIM_simInitList$noPM$modelBiomass$mod

  pars <- unlist(getME(modBiomass, c("theta")))
  updateModBiomass <- update(modBiomass, devFunOnly = TRUE, data = modBiomass@frame)
  grad <- numDeriv::grad(updateModBiomass, pars)
  hess <- numDeriv::hessian(updateModBiomass, pars)
  sc_grad <- solve(hess, grad)

  if (length(modBiomass@optinfo$conv$lme4$messages) &
      max(pmin(abs(sc_grad), abs(grad))) > 0.001) {
    ## Nelder and L-BFGS-B methods didn't converge
    # , REML = FALSE,
    # control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb')))
    # ss <- getME(modBiomass, c("theta","fixef"))
    # modBiomass <- update(modBiomass, start = ss, data = modBiomass@frame,
    #                            control = lmerControl(calc.derivs = FALSE,  ## it's faster
    #                                                  optimizer = "optimx",
    #                                                  optCtrl = list(method='nlminb', maxit = 1e5)))  ## keeps hitting maxit...
    ncores <- detectCores()
    library(dfoptim)
    .specialData <<- modBiomass@frame
    diff_optims <- allFit(modBiomass$mod, parallel = 'multicore', ncpus = ncores)
    is.OK <- sapply(diff_optims, is, "merMod")  ## nlopt NELDERMEAD failed, others succeeded
    aa.OK <- diff_optims[is.OK]
    lapply(aa.OK,function(x) x@optinfo$conv$lme4$messages)

    pars <- unlist(getME(modBiomass, c("theta")))
    grad <- numDeriv::grad(update(modBiomass, devFunOnly = TRUE), pars)
    hess <- numDeriv::hessian(update(modBiomass, devFunOnly = TRUE), pars)
    sc_grad <- solve(hess, grad)
    max(pmin(abs(sc_grad),abs(grad)))
  }
}


## DIAGNOSE FIRE IGNITION MODEL ----------------------------------
## to get predicted values we re-run the IgnitionPredict with the original data
if (FALSE) {
  library(data.table)
  library(ggplot2)
  library(ggspatial)
  library(ggpubr)

  objects <- list(
    "fireSense_IgnitionFitted" = LIM_simInitList$noPM$fireSense_IgnitionFitted
    , "dataFireSense_IgnitionFit" = LIM_simInitList$noPM$dataFireSense_IgnitionFit
  )

  parameters <- list(
    fireSense_IgnitionPredict = list(
      "data" = "dataFireSense_IgnitionFit"
      , "modelObjName" = "fireSense_IgnitionFitted" # This is the default
      , "rescaleFactor" = 1
    )
  )
  simOutFireFreqPredVals <- Cache(simInitAndSpades
                                  , times = list(start = 0, end = 0)
                                  , modules = "fireSense_IgnitionPredict"
                                  , paths = simPaths
                                  , objects = objects
                                  , params = parameters
                                  , cacheRepo = simPaths$cachePath
                                  , userTags = "fireSense_IgnitionPredict"
                                  , omitArgs = c("userTags")
  )

  ## plot predicted and fitted values
  ignitionsData <- LIM_simInitList$noPM$dataFireSense_IgnitionFit
  ignitionsData[,  rows := 1:nrow(ignitionsData)]
  predVals <- data.table(rows = as.integer(names(simOutFireFreqPredVals$fireSense_IgnitionPredicted)),
                         fittedVals = simOutFireFreqPredVals$fireSense_IgnitionPredicted)
  ignitionsData <- predVals[ignitionsData, on = .(rows)]
  ignitionsData[, n_firesPred := rpois(.N, fittedVals)]

  plot1 <- ggplot(data = ignitionsData, aes(y = fittedVals, x = n_fires)) +
    geom_point() +
    labs(y = "lambda", x = "observed no. fires")

  plot2 <- ggplot(data = ignitionsData, aes(y = fittedVals, x = n_firesPred)) +
    geom_point() +
    labs(y = "lambda", x = "predicted (rpois) no. fires")

  plot3 <- ggplot(data = ignitionsData, aes(y = n_firesPred, x = n_fires)) +
    geom_point() +
    labs(y = "predicted (rpois) no. fires", x = "observed no. fires")

  # plotData <- ignitionsData[, list(obsFires = sum(n_fires), predFires = sum(n_firesPred)),
  #                           by = year]
  plotData <- ignitionsData[, list(obsFires = sum(n_fires, na.rm = TRUE),
                                   predFires = sum(fittedVals, na.rm = TRUE)),
                            by = year]
  plotData <- melt(plotData, id.var = "year")
  plot4 <- ggplot(data = plotData, aes(y = value, x = year, colour = variable)) +
    geom_line(size = 1) +
    scale_color_discrete(labels = c("obsFires" = "sum of observed no. fires",
                                    "predFires" = "sum of predicted (lambda) no. fires")) +
    theme_pubr(margin = FALSE, base_size = 14) +
    theme(legend.position = "bottom") +
    labs(y = "no. fires", x = "year", colour = "")

  gridExtra::grid.arrange(plot1, plot2, plot3, plot4)

  plot5 <- ggplot() +
    layer_spatial(data = simOutFireFreq$fireSense_IgnitionPredicted) +
    layer_spatial(data = simOutPreSim$fireLocations, colour = "darkred") +
    annotation_north_arrow(style = north_arrow_minimal,
                           location = "tr", which_north = "true") +
    scale_fill_distiller(palette = "Greys", na.value = "transparent", direction = 1) +
    theme_pubr(margin = FALSE, legend = "right", base_size = 14) +
    labs(x = "longitude", y = "latitude", fill = expression(lambda))

  plot6 <- ggarrange(plot4, plot5, widths = c(0.6,0.4),
                     labels = "auto", font.label = list(size = 20))
  ggsave("C:/Users/Ceres Barros/Google Drive/Shared/McIntire-lab/Manuscripts_inPrep/LIMmodel_paper",
         plot6, width = 12, height = 7, dpi = 300)
}

