## this has been moved to ToolsCB


# ## ------------------------------------
# ## LATEX/RMARKDOWN CORRELATION TABLES
# ## ------------------------------------
#
# ## function to make a nice correlation table for latex/Rmd
# ## usage example:
# # xtable(prepCorrTable(DT))
#
# ## adapted from http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html
# require(Hmisc)
#
# ## x matrix or any object compatible with as.matrix, from where the covariance matrix will be calculated
# ## method see Hmisc::rcorr 'type' argument.
# prepCorrTable <- function(x, method = "pearson") {
#   ## checks
#   if (!all(class(x) == "matrix"))
#     x <- as.matrix(x)
#   if (!method %in% c("pearson", "spearman"))
#     stop("method must be 'pearson' or 'spearman'")
#
#   corValue <- rcorr(x, type = method)$r
#   p <- rcorr(x, type = method)$P
#
#   ## define notions for significance levels; spacing is important.
#   mystars <- ifelse(p < .001, "***", ifelse(p < .01, "** ", ifelse(p < .05, "* ", " ")))
#
#   ## trunctuate the matrix that holds the correlations to two decimal
#   corValue <- format(round(cbind(rep(-1.11, ncol(x)), corValue), 2))[,-1]
#
#   ## build a new matrix that includes the correlations with their apropriate stars
#   Rnew <- matrix(paste(corValue, mystars, sep=""), ncol=ncol(x))
#   diag(Rnew) <- paste(diag(corValue), " ", sep="")
#   rownames(Rnew) <- colnames(x)
#   colnames(Rnew) <- paste(colnames(x), "", sep="")
#
#   ## remove upper triangle
#   Rnew <- as.matrix(Rnew)
#   Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
#   Rnew <- as.data.frame(Rnew)
#
#   ## remove last column and return the matrix (which is now a data frame)
#   Rnew <- cbind(Rnew[1:length(Rnew)-1])
#   return(Rnew)
# }
