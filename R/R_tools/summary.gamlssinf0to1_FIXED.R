## --------------------------------------
## SUMMARY FOR GAMLSSINF0TO1 FUNCTION
## --------------------------------------

## this is a fix for the gamlss.inf:::summary.gamlssinf0to1 function which has issues when
## random terms were included in the model.
## because these random terms have NA's in the covariante matrix, the t-values and p-values
## would not be calcualted, and the summary would be screwed up (like variables were missing,
## when the random terms were simply missing, but names were not updated)


summary.gamlssinf0to1_2 <- function (object, type = c("vcov", "qr"), robust = FALSE, save = FALSE,
          hessian.fun = c("R", "PB"), digits = max(3, getOption("digits") -
                                                     3), ...) {
  type <- match.arg(type)
  pm <- ps <- pn <- pt <- px0 <- px1 <- 0
  mu.coef.table <- NULL
  sigma.coef.table <- NULL
  nu.coef.table <- NULL
  tau.coef.table <- NULL
  xi0.coef.table <- NULL
  xi1.coef.table <- NULL
  if (type == "vcov") {
    covmat <- try(suppressWarnings(vcov.gamlssinf0to1(object,
                                                      type = "all", robust = robust, hessian.fun = hessian.fun)),
                  silent = TRUE)
    if (any(class(covmat) %in% "try-error" || any(is.na(covmat$se)))) {
      warning(paste("summary: vcov has failed, option qr is used instead\n"))
      type <- "qr"
    }
  }
  ifWarning <- rep(FALSE, length(object$parameters))
  if (type == "vcov") {
    ## exclude random effects
    coef <- covmat$coef[!grepl("random\\(|re\\(", names(covmat$coef))]
    se <- covmat$se[!grepl("random\\(|re\\(", names(covmat$se))]
    tvalue <- coef/se
    pvalue <- 2 * pt(-abs(tvalue), object$df.res)
    coef.table <- cbind(coef, se, tvalue, pvalue)
    dimnames(coef.table) <- list(names(coef), c("Estimate",
                                                "Std. Error", "t value", "Pr(>|t|)"))
    cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n")
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n")
    est.disp <- FALSE
    if ("mu" %in% object$parameters) {
      ifWarning[1] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "mu"), specials = .gamlss.sm.list), "specials"))))
      if (object$mu.df != 0) {
        pm <- object$mu.qr$rank
        p1 <- 1:pm
        cat("-------------------------------------------------------------------\n")
        cat("Mu link function: ", object$mu.link)
        cat("\n")
        cat("Mu Coefficients:")
        if (is.character(co <- object$contrasts))
          cat("  [contrasts: ", apply(cbind(names(co),
                                            co), 1, paste, collapse = "="), "]")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$mu.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Mu parameter is fixed \n")
        if (all(object$mu.fv == object$mu.fv[1]))
          cat("Mu = ", object$mu.fv[1], "\n")
        else cat("Mu is equal with the vector (", object$mu.fv[1],
                 ",", object$mu.fv[2], ",", object$mu.fv[3],
                 ",", object$mu.fv[4], ", ...) \n")
      }
    }
    if ("sigma" %in% object$parameters) {
      ifWarning[2] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "sigma"), specials = .gamlss.sm.list), "specials"))))
      if (object$sigma.df != 0) {
        ps <- object$sigma.qr$rank
        p1 <- (pm + 1):(pm + ps)
        cat("-------------------------------------------------------------------\n")
        cat("Sigma link function: ", object$sigma.link)
        cat("\n")
        cat("Sigma Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$sigma.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Sigma parameter is fixed")
        cat("\n")
        if (all(object$sigma.fv == object$sigma.fv[1]))
          cat("Sigma = ", object$sigma.fv[1], "\n")
        else cat("Sigma is equal with the vector (",
                 object$sigma.fv[1], ",", object$sigma.fv[2],
                 ",", object$sigma.fv[3], ",", object$sigma.fv[4],
                 ", ...) \n")
      }
    }
    if ("nu" %in% object$parameters) {
      ifWarning[3] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "nu"), specials = .gamlss.sm.list), "specials"))))
      if (object$nu.df != 0) {
        pn <- object$nu.qr$rank
        p1 <- (pm + ps + 1):(pm + ps + pn)
        cat("-------------------------------------------------------------------\n")
        cat("Nu link function: ", object$nu.link, "\n")
        cat("Nu Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$nu.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Nu parameter is fixed")
        cat("\n")
        if (all(object$nu.fv == object$nu.fv[1]))
          cat("Nu = ", object$nu.fv[1], "\n")
        else cat("Nu is equal with the vector (", object$nu.fv[1],
                 ",", object$nu.fv[2], ",", object$nu.fv[3],
                 ",", object$nu.fv[4], ", ...) \n")
      }
    }
    if ("tau" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "tau"), specials = .gamlss.sm.list), "specials"))))
      if (object$tau.df != 0) {
        pt <- object$tau.qr$rank
        p1 <- (pm + ps + pn + 1):(pm + ps + pn + pt)
        cat("-------------------------------------------------------------------\n")
        cat("Tau link function: ", object$tau.link,
            "\n")
        cat("Tau Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$tau.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Tau parameter is fixed")
        cat("\n")
        if (all(object$tau.fv == object$tau.fv[1]))
          cat("Tau = ", object$tau.fv[1], "\n")
        else cat("Tau is equal with the vector (", object$tau.fv[1],
                 ",", object$tau.fv[2], ",", object$tau.fv[3],
                 ",", object$tau.fv[4], ", ...) \n")
      }
    }
    if ("xi0" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "xi0"), specials = .gamlss.sm.list), "specials"))))
      if (object$xi0.df != 0) {
        px0 <- object$xi0.qr$rank
        p1 <- (pm + ps + pn + pt + 1):(pm + ps + pn +
                                         pt + px0)
        cat("-------------------------------------------------------------------\n")
        cat("xi0 link function: ", object$xi0.link,
            "\n")
        cat("xi0 Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$xi0.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("xi0 parameter is fixed")
        cat("\n")
        if (all(object$xi0.fv == object$xi0.fv[1]))
          cat("xi0 = ", object$xi0.fv[1], "\n")
        else cat("xi0 is equal with the vector (", object$xi0.fv[1],
                 ",", object$xi0.fv[2], ",", object$xi0.fv[3],
                 ",", object$xi0.fv[4], ", ...) \n")
      }
    }
    else px0 <- 0
    if ("xi1" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "xi1"), specials = .gamlss.sm.list), "specials"))))
      if (object$xi1.df != 0) {
        px1 <- object$xi1.qr$rank
        p1 <- (pm + ps + pn + pt + px0 + 1):(pm + ps +
                                               pn + pt + px0 + px1)
        cat("-------------------------------------------------------------------\n")
        cat("xi1 link function: ", object$xi1.link,
            "\n")
        cat("xi1 Coefficients:")
        cat("\n")
        printCoefmat(coef.table[p1, , drop = FALSE],
                     digits = digits, signif.stars = TRUE)
        cat("\n")
      }
      else if (object$xi1.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("xi1 parameter is fixed")
        cat("\n")
        if (all(object$xi1.fv == object$xi1.fv[1]))
          cat("xi1 = ", object$xi1.fv[1], "\n")
        else cat("xi1 is equal with the vector (", object$xi1.fv[1],
                 ",", object$xi1.fv[2], ",", object$xi1.fv[3],
                 ",", object$xi1.fv[4], ", ...) \n")
      }
    }
    if (any(ifWarning)) {
      cat("-------------------------------------------------------------------\n")
      cat("NOTE: Additive smoothing terms exist in the formulas: \n")
      cat(" i) Std. Error for smoothers are for the linear effect only. \n")
      cat("ii) Std. Error for the linear terms maybe are not accurate. \n")
    }
    cat("-------------------------------------------------------------------\n")
    cat("No. of observations in the fit: ", object$noObs,
        "\n")
    cat("Degrees of Freedom for the fit: ", object$df.fit)
    cat("\n")
    cat("      Residual Deg. of Freedom: ", object$df.residual,
        "\n")
    cat("                      at cycle: ", object$iter,
        "\n \n")
    cat("Global Deviance:    ", object$G.deviance, "\n            AIC:    ",
        object$aic, "\n            SBC:    ", object$sbc,
        "\n")
    cat("*******************************************************************")
    cat("\n")
  }
  if (type == "qr") {
    estimatesgamlss <- function(object, Qr, p1, coef.p,
                                est.disp, df.r, digits = max(3, getOption("digits") -
                                                               3), covmat.unscaled, ...) {
      dimnames(covmat.unscaled) <- list(names(coef.p),
                                        names(coef.p))
      covmat <- covmat.unscaled
      var.cf <- diag(covmat)
      s.err <- sqrt(var.cf)
      tvalue <- coef.p/s.err
      dn <- c("Estimate", "Std. Error")
      if (!est.disp) {
        pvalue <- 2 * pnorm(-abs(tvalue))
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p),
                                     c(dn, "z value", "Pr(>|z|)"))
      }
      else if (df.r > 0) {
        pvalue <- 2 * pt(-abs(tvalue), df.r)
        coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
        dimnames(coef.table) <- list(names(coef.p),
                                     c(dn, "t value", "Pr(>|t|)"))
      }
      else {
        coef.table <- cbind(coef.p, Inf)
        dimnames(coef.table) <- list(names(coef.p),
                                     dn)
      }
      return(coef.table)
    }
    dispersion <- NULL
    cat("*******************************************************************")
    cat("\nFamily: ", deparse(object$family), "\n")
    cat("\nCall: ", deparse(object$call), "\n", fill = TRUE)
    cat("Fitting method:", deparse(object$method), "\n\n")
    est.disp <- FALSE
    df.r <- object$noObs - object$mu.df
    if ("mu" %in% object$parameters) {
      ifWarning[1] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "mu"), specials = .gamlss.sm.list), "specials"))))
      if (object$mu.df != 0) {
        Qr <- object$mu.qr
        df.r <- object$noObs - object$mu.df
        if (is.null(dispersion))
          dispersion <- if (any(object$family == c("PO",
                                                   "BI", "EX", "P1")))
            1
        else if (df.r > 0) {
          est.disp <- TRUE
          if (any(object$weights == 0))
            warning(paste("observations with zero weight",
                          "not used for calculating dispersion"))
        }
        else Inf
        p <- object$mu.df
        p1 <- 1:(p - object$mu.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        mu.coef.table <- estimatesgamlss(object = object,
                                         Qr = object$mu.qr, p1 = p1, coef.p = object$mu.coefficients[Qr$pivot[p1]],
                                         est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("Mu link function: ", object$mu.link)
        cat("\n")
        cat("Mu Coefficients:")
        if (is.character(co <- object$contrasts))
          cat("  [contrasts: ", apply(cbind(names(co),
                                            co), 1, paste, collapse = "="), "]")
        cat("\n")
        printCoefmat(mu.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$mu.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Mu parameter is fixed")
        cat("\n")
        if (all(object$mu.fv == object$mu.fv[1]))
          cat("Mu = ", object$mu.fv[1], "\n")
        else cat("Mu is equal with the vector (", object$mu.fv[1],
                 ",", object$mu.fv[2], ",", object$mu.fv[3],
                 ",", object$mu.fv[4], ", ...) \n")
      }
      coef.table <- mu.coef.table
    }
    else {
      if (df.r > 0) {
        est.disp <- TRUE
        if (any(object$weights == 0))
          warning(paste("observations with zero weight",
                        "not used for calculating dispersion"))
      }
    }
    if ("sigma" %in% object$parameters) {
      ifWarning[2] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "sigma"), specials = .gamlss.sm.list), "specials"))))
      if (object$sigma.df != 0) {
        Qr <- object$sigma.qr
        df.r <- object$noObs - object$sigma.df
        p <- object$sigma.df
        p1 <- 1:(p - object$sigma.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        sigma.coef.table <- estimatesgamlss(object = object,
                                            Qr = object$sigma.qr, p1 = p1, coef.p = object$sigma.coefficients[Qr$pivot[p1]],
                                            est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("Sigma link function: ", object$sigma.link)
        cat("\n")
        cat("Sigma Coefficients:")
        cat("\n")
        printCoefmat(sigma.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$sigma.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Sigma parameter is fixed")
        cat("\n")
        if (all(object$sigma.fv == object$sigma.fv[1]))
          cat("Sigma = ", object$sigma.fv[1], "\n")
        else cat("Sigma is equal with the vector (",
                 object$sigma.fv[1], ",", object$sigma.fv[2],
                 ",", object$sigma.fv[3], ",", object$sigma.fv[4],
                 ", ...) \n")
      }
      coef.table <- rbind(mu.coef.table, sigma.coef.table)
    }
    if ("nu" %in% object$parameters) {
      ifWarning[3] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "nu"), specials = .gamlss.sm.list), "specials"))))
      if (object$nu.df != 0) {
        Qr <- object$nu.qr
        df.r <- object$noObs - object$nu.df
        p <- object$nu.df
        p1 <- 1:(p - object$nu.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        nu.coef.table <- estimatesgamlss(object = object,
                                         Qr = object$nu.qr, p1 = p1, coef.p = object$nu.coefficients[Qr$pivot[p1]],
                                         est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("Nu link function: ", object$nu.link, "\n")
        cat("Nu Coefficients:")
        cat("\n")
        printCoefmat(nu.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$nu.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Nu parameter is fixed")
        cat("\n")
        if (all(object$nu.fv == object$nu.fv[1]))
          cat("Nu = ", object$nu.fv[1], "\n")
        else cat("Nu is equal with the vector (", object$nu.fv[1],
                 ",", object$nu.fv[2], ",", object$nu.fv[3],
                 ",", object$nu.fv[4], ", ...) \n")
      }
      coef.table <- rbind(mu.coef.table, sigma.coef.table,
                          nu.coef.table)
    }
    if ("tau" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "tau"), specials = .gamlss.sm.list), "specials"))))
      if (object$tau.df != 0) {
        Qr <- object$tau.qr
        df.r <- object$noObs - object$tau.df
        p <- object$tau.df
        p1 <- 1:(p - object$tau.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        tau.coef.table <- estimatesgamlss(object = object,
                                          Qr = object$tau.qr, p1 = p1, coef.p = object$tau.coefficients[Qr$pivot[p1]],
                                          est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("Tau link function: ", object$tau.link,
            "\n")
        cat("Tau Coefficients:")
        cat("\n")
        printCoefmat(tau.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$tau.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("Tau parameter is fixed")
        cat("\n")
        if (all(object$tau.fv == object$tau.fv[1]))
          cat("Tau = ", object$tau.fv[1], "\n")
        else cat("Tau is equal with the vector (", object$tau.fv[1],
                 ",", object$tau.fv[2], ",", object$tau.fv[3],
                 ",", object$tau.fv[4], ", ...) \n")
      }
      coef.table <- rbind(mu.coef.table, sigma.coef.table,
                          nu.coef.table, tau.coef.table)
    }
    if ("xi0" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "xi0"), specials = .gamlss.sm.list), "specials"))))
      if (object$xi0.df != 0) {
        Qr <- object$xi0.qr
        df.r <- object$noObs - object$xi0.df
        p <- object$xi0.df
        p1 <- 1:(p - object$xi0.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        xi0.coef.table <- estimatesgamlss(object = object,
                                          Qr = object$xi0.qr, p1 = p1, coef.p = object$xi0.coefficients[Qr$pivot[p1]],
                                          est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("xi0 link function: ", object$xi0.link,
            "\n")
        cat("xi0 Coefficients:")
        cat("\n")
        printCoefmat(xi0.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$xi0.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("xi0 parameter is fixed")
        cat("\n")
        if (all(object$xi0.fv == object$xi0.fv[1]))
          cat("xi0 = ", object$xi0.fv[1], "\n")
        else cat("xi0 is equal with the vector (", object$xi0.fv[1],
                 ",", object$xi0.fv[2], ",", object$xi0.fv[3],
                 ",", object$xi0.fv[4], ", ...) \n")
      }
      coef.table <- rbind(mu.coef.table, sigma.coef.table,
                          nu.coef.table, tau.coef.table, xi0.coef.table)
    }
    if (!("xi0" %in% object$parameters))
      xi0.coef.table <- NULL
    if ("xi1" %in% object$parameters) {
      ifWarning[4] <- (!is.null(unlist(attr(terms(formula(object,
                                                          "xi1"), specials = .gamlss.sm.list), "specials"))))
      if (object$xi1.df != 0) {
        Qr <- object$xi1.qr
        df.r <- object$noObs - object$xi1.df
        p <- object$xi1.df
        p1 <- 1:(p - object$xi1.nl.df)
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        xi1.coef.table <- estimatesgamlss(object = object,
                                          Qr = object$xi1.qr, p1 = p1, coef.p = object$xi1.coefficients[Qr$pivot[p1]],
                                          est.disp = est.disp, df.r = df.r, covmat.unscaled = covmat.unscaled)
        cat("-------------------------------------------------------------------\n")
        cat("xi1 link function: ", object$xi1.link,
            "\n")
        cat("xi1 Coefficients:")
        cat("\n")
        printCoefmat(xi1.coef.table, digits = digits,
                     signif.stars = TRUE)
        cat("\n")
      }
      else if (object$xi1.fix == TRUE) {
        cat("-------------------------------------------------------------------\n")
        cat("xi1 parameter is fixed")
        cat("\n")
        if (all(object$xi1.fv == object$xi1.fv[1]))
          cat("xi1 = ", object$xi1.fv[1], "\n")
        else cat("xi1 is equal with the vector (", object$xi1.fv[1],
                 ",", object$xi1.fv[2], ",", object$xi1.fv[3],
                 ",", object$xi1.fv[4], ", ...) \n")
      }
      coef.table <- rbind(mu.coef.table, sigma.coef.table,
                          nu.coef.table, tau.coef.table, xi0.coef.table,
                          xi1.coef.table)
    }
    if (!("xi0" %in% object$parameters))
      xi0.coef.table <- NULL
    if (any(ifWarning)) {
      cat("-------------------------------------------------------------------\n")
      cat("NOTE: Additive smoothing terms exist in the formulas: \n")
      cat(" i) Std. Error for smoothers are for the linear effect only. \n")
      cat("ii) Std. Error for the linear terms may not be reliable. \n")
    }
    cat("-------------------------------------------------------------------\n")
    cat("No. of observations in the fit: ", object$noObs,
        "\n")
    cat("Degrees of Freedom for the fit: ", object$df.fit)
    cat("\n")
    cat("      Residual Deg. of Freedom: ", object$df.residual,
        "\n")
    cat("                      at cycle: ", object$iter,
        "\n \n")
    cat("Global Deviance:    ", object$G.deviance, "\n            AIC:    ",
        object$aic, "\n            SBC:    ", object$sbc,
        "\n")
    cat("*******************************************************************")
    cat("\n")
  }
  if (save == TRUE) {
    out <- as.list(environment())
    return(out)
  }
  invisible(coef.table)
}
