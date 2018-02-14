
# Everything in this file gets sourced during simInit, and all functions and objects
# are put into the simList. To use objects and functions, use sim$xxx.
defineModule(sim, list(
  name = "fireStats",
  description = "Basic statistical analysis/visuals from fire state-and-transition [toy] model results",
  keywords = c("fire", "basic statistics"),
  authors = person("Ceres", "Barros", email = "cbarros@mail.ubc.ca", role = c("aut", "cre")),
  childModules = character(0),
  version = numeric_version("0.0.1"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "fire_STSM.Rmd"),
  reqdPkgs = list("raster", "data.table", "ggplot2", "reshape2"),
  parameters = rbind(
    defineParameter(".plotStats", "logical", TRUE, NA, NA, "This describes whether sumamry statistics should be plotted or not")
  ),
  inputObjects = bind_rows(
    expectsInput(objectName = NA, objectClass = NA, desc = NA, sourceURL = NA)
  ),
  outputObjects = bind_rows(
    createsOutput(objectName = NA, objectClass = NA, desc = NA)
  )
))


doEvent.fireStats = function(sim, eventTime, eventType, debug = FALSE) {
  switch(
    eventType,
    init = {
      ## do stuff for this event
      sim <- statsInit(sim)
      
      ## schedule plotting event
      if(P(sim)$.plotStats) {
        sim <- scheduleEvent(sim, eventTime = end(sim), moduleName = "fireStats", eventType = "dataPrep", eventPriority = 5)
        sim <- scheduleEvent(sim, eventTime = end(sim), moduleName = "fireStats", eventType = "plot", eventPriority = 6)
      }
    },
    dataPrep = {
      ## do stuff for this event
      sim <- statsPrep(sim)
      
    },
    plot = {
      ## do stuff for this event
      sim <- statsPlot(sim)

    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
  return(invisible(sim))
}


## Module initialization - nothing to be done, so left empty
statsInit <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #

  # ! ----- STOP EDITING ----- ! #

  return(invisible(sim))
}

## Plots of basic statistics 
statsPrep <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  
  ## Calculate  relative abundances of vegetation types per year
  freqs <- lapply(sim$vegetation, FUN = function(x) {
    return(table(x[]))
  })
  
  relativeAbundance <- do.call(cbind, lapply(freqs, FUN = function(x) {
    x = as.matrix(x)
    relativeAbundance <- apply(x, 1, FUN = function(y) y/colSums(x))
    return(as.matrix(relativeAbundance))
  }))

  relativeAbundance <- melt(relativeAbundance)
  names(relativeAbundance) = c("VegType", "Year", "Abundance")
  relativeAbundance$VegType <- factor(relativeAbundance$VegType, levels = paste0(0:5),
                                      labels = c("Non-burnable", "Grass", "Shrub", "Deciduous", "Mixed", "Coniferous"))
  
  sim$vegRelativAbunds <- relativeAbundance
  
  ## Calculate average fire size per year
  sim$fireSizes <- as.data.frame(sapply(sim$spreadRas, FUN = function(ras) {
    fires <- na.exclude(unique(ras[]))
    fire_sizes <- sapply(fires, FUN = function(x) {
      sum(ras[]==x, na.rm = TRUE)
    })
    return(fire_sizes)
  }))
  
  sim$fireSizes$Fire_ID <- rownames(sim$fireSizes)
  sim$fireSizes <- melt(sim$fireSizes, id.vars = "Fire_ID", value.name = "Size", variable.name = "Year")
  sim$fireSizes$Year <- as.numeric(as.character(sub("V", "", sim$fireSizes$Year)))
  
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}

## Plots of basic statistics 
statsPlot <- function(sim) {
  # ! ----- EDIT BELOW ----- ! #
  print(
    ggplot() +
      geom_line(data = sim$vegRelativAbunds, aes(x = Year, y = Abundance, col = VegType), size = 1) +
      stat_summary(data = sim$fireSizes, aes(x = Year, y = Size/1000), fun.y = "mean", 
                   col = "red", geom = "line", size = 1) +
      stat_summary(data = sim$fireSizes, aes(x = Year, y = Size/1000), fun.data = "mean_se", 
                   col = "red") +
      theme_bw() +
      labs(title = "Changes in vegetation\nfire size (/100) in red")
  )
  # ! ----- STOP EDITING ----- ! #
  return(invisible(sim))
}



.inputObjects = function(sim) {
  
  return(invisible(sim))
}
### add additional events as needed by copy/pasting from above
