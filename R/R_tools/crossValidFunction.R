## ----------------------------------------
## RUN K-FOLD CROSS VALIDATION
##
## Ceres Feb 26th 2020
## ----------------------------------------

## this script should be sourced

## CROSS-CALIDATION FUNCTION

## fullDT - data.table with full dataset
## statsModel - the statistical model to validate. Only works with gamlss models
## k - integer with number of chunks that the data should be partioned in
## idCol - column with pixel/observation IDs (optional)
## origDataVars - the data used to fit the statsModule, needs to be passed to gamlss::predictAll
##   (it may not be able to access it) bu also to make sure newdata in gamlss::predictAll
##  has the same variables (even if they're not used in the model)
## cacheObj1/2 - an object used by Cache for digesting, to avoid digesting the (potentially) large data arguments

crossValidFunction <- function (fullDT, statsModel, origData, k = 4, idCol, cacheObj1, cacheObj2) {
  if (!is.null(idCol))
    origDataVars <- c(names(origData), idCol)

  ## remove NAs from the data without subsetting columns
  toKeep <- na.omit(fullDT[, ..origDataVars])[, ..idCol]
  setkeyv(fullDT, idCol)
  setkeyv(toKeep, idCol)
  fullDT <- fullDT[toKeep]

  ## make chunks of 1/4 of the data
  cols2 <- c("FIRE_NAME", idCol)
  sampDT <- fullDT[,..cols2]
  set.seed(123)
  sampDT[, sampID := sample(1:k, size = length(get(idCol)), replace = TRUE),
         by = FIRE_NAME]
  ## join samp IDs with data
  fullDT <- sampDT[fullDT, on = cols2]
  rm(cols2)

  origDataVars <- c(origDataVars, "sampID")

  crossValidResults <- lapply(unique(fullDT$sampID), FUN = calcCrossValidMetrics,
                              fullDT = fullDT, origData = origData, idCol = idCol,
                              statsModel = statsModel, origDataVars = origDataVars)

  return(crossValidResults)
}


## CALCULATE VALIDATION METRICS AND CONFUSION MATRIX
## to allow caching without digesting the large data table

## samp - the sample number to pick to use as the test data set
## fullDT - the full dataset (not necessarily the one used to fit `statsModel`, which could have been a subset (e.g. fewer columns))
## origData - the data used to fit `statsModel`
## statsModel - the fitted model
## origDataVars - a character vector of the variables used in model fitting (including response variable and random effects.)


## outputs a list with 2 entries
calcCrossValidMetrics <- function(samp, fullDT, origData, statsModel, origDataVars) {
  ## predict requires the original and new data to have the same columns
  if (!all(names(origData) %in% names(fullDT)))
    stop("'fullDT' needs to include all the columns in 'origData'")

  ## subset
  trainData <- na.omit(fullDT[sampID != samp, ..origDataVars])
  testData <- na.omit(fullDT[sampID == samp, ..origDataVars])

  ## checks
  if (length(setdiff(unique(fullDT$FIRE_NAME),
                     unique(testData$FIRE_NAME))) |
      length(setdiff(unique(fullDT$FIRE_NAME),
                     unique(trainData$FIRE_NAME))))
    stop("Fires lost in sampling!")


  ## trainData an testData cannot have extra cols
  cols <- names(origData)
  trainData <<- trainData[, ..cols]
  testData <- testData[, ..cols]

  ## refit model on training sample - this is failing due to singularity
  ## then predict
  trainModel <- update(statsModel, data = trainData)
  predictionsDT <- predictAll(trainModel,
                              newdata = testData, data = trainData,
                              type = "response", output = "matrix")
  predictionsDT <- as.data.table(predictionsDT)

  ## change name
  setnames(predictionsDT, "y", "SEV_PROP")

  ## predict using rBEINF
  ## tried generating many values and averaging, but doing that results in the same value
  predictionsDT[, predSEV_PROP := mean(rBEINF(10, mu, sigma, nu, tau)),
                by = row.names(predictionsDT)]

  ## add severity classes
  testData <- na.omit(fullDT[sampID == samp, ..origDataVars]) ## redo testData in case idCol was dropped when subsetting to model data
  predictionsDT[, pixID := testData$pixID]
  predictionsDT <- fullDT[, .(pixID, SEV_CLASS)][predictionsDT, on = "pixID"]

  ## convert to classes, using the quantiles corresponding to the observed class proportions
  ## accumulate proportions to get probabilities
  quantProbs <- cumsum(table(predictionsDT$SEV_CLASS)/nrow(predictionsDT))
  classRanges <- c(0, quantile(predictionsDT$predSEV_PROP, probs = quantProbs))

  predictionsDT[, predSEV_CLASS := cut(predSEV_PROP, breaks = classRanges,
                                       include.lowest = TRUE, right = FALSE)]  ## classify as with intervals as ],]

  ## convert to numbered factor (subtracting one, because classes are 0-5)
  predictionsDT[, predSEV_CLASS := as.numeric(predSEV_CLASS)-1]
  classes <- as.character(sort(unique(predictionsDT$SEV_CLASS)))
  predictionsDT[, `:=`(SEV_CLASS = factor(SEV_CLASS, levels = classes),
                       predSEV_CLASS = factor(predSEV_CLASS, levels = classes))]

  ## VALIDATION STATISTICS WITH CLASSES ----------------------------------
  ## calculate overall statistics
  validMetrics <- caret::multiClassSummary(predictionsDT[, list(obs = SEV_CLASS,
                                                                pred = predSEV_CLASS)],
                                           lev = classes)
  ## calculate confusion matrix
  confMatrix <- caret::confusionMatrix(data = predictionsDT$predSEV_CLASS, reference = predictionsDT$SEV_CLASS)

  ## VALIDATION STATISTICS WITH CONTINUOUS VARIABLE -----------------------
  Rsquared <- caret::postResample(pred = predictionsDT$predSEV_PROP, obs = predictionsDT$SEV_PROP)
  Rsquared <- Rsquared["Rsquared"]

  ## COEFFICIENTS
  list(validMetrics = validMetrics, confMatrix = confMatrix,
       Rsquared = Rsquared, coefs = coefAll(trainModel))
}
