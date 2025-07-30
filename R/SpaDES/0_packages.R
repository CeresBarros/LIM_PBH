## ------------------------------------------------------
## FIRE MODELLING WITH SpaDES -- package installation never used
##
## Ceres: July 2025
## ------------------------------------------------------

## note that on coco Docker I had to:
install.packages("http://cran.r-project.org/src/contrib/Archive/sp/sp_1.4-6.tar.gz")
install.packages("http://cran.r-project.org/src/contrib/Archive/rgdal/rgdal_1.5-30.tar.gz")
install.packages("http://cran.r-project.org/src/contrib/Archive/terra/terra_1.5-34.tar.gz")
install.packages("http://cran.r-project.org/src/contrib/Archive/sf/sf_1.0-7.tar.gz")
devtools::install_github("CeresBarros/ToolsCB@4b12ff37e29455350abbaa340f1fedc6aee38a49", dependencies = FALSE)


## CODE BELOW NEVER USED ON DOCKER:
## before switching to proj lib, install an updated version of Require in default lib
remotes::install_github("PredictiveEcology/Require@70720ac5fd104b37c0c72f7613d79b81ec437dc7")

## RESTART R AND KEEP GOING

library(Require); library(httr2) ## load from default/non-proj library

if (!exists("pkgDir")) {
  pkgDir <- file.path(
    if (Sys.info()[["user"]] %in% c("rstudio", "root")) "packages_docker" else "packages",
    version$platform,
    paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
  )

  if (!dir.exists(pkgDir)) {
    dir.create(pkgDir, recursive = TRUE)
  }
  .libPaths(pkgDir)
}

## old code
# if (!require("Require")) {
#   # remotes::install_github("PredictiveEcology/Require@7eaa3af6443fa9acf8ef461d9a02e544174eda38", upgrade = FALSE)
# }

if (FALSE) {
  ## write file
  # Require::pkgSnapshot("packages_docker/pkgSnapshot.txt", standAlone = TRUE, exact = TRUE)
  # Require::pkgSnapshot("packages/pkgSnapshot.txt", standAlone = TRUE, exact = TRUE)

  # Much later on a different or same machine
  # Require::Require(packageVersionFile = "packages_docker/pkgSnapshot.txt", standAlone = TRUE)
  Require(packageVersionFile = "packages/pkgSnapshot.txt", standAlone = TRUE)

  ## these versions need to be ensured and if not on the snapshot file
  Install("PredictiveEcology/Require@7eaa3af6443fa9acf8ef461d9a02e544174eda38", dependencies = FALSE)
  Install("PredictiveEcology/reproducible@97f147033f10061dfe40da4ce76575dfef89dfbd", dependencies = FALSE)
  Install("CeresBarros/SpaDES.core@2e8b95b04cf6b93bc45796ddc8a55c8d1618432d", dependencies = FALSE)
  Install("CeresBarros/SpaDES.tools@780cd50cbf156faa86d7cb3a609c6d5edf359752", dependencies = FALSE)
  Install("ianmseddy/PSPclean@3c3f0e7082e14c111a607c3ba803abf0396343e6", dependencies = FALSE)
  Install("CeresBarros/SpaDES.experiment@82ffdb7cce912013e19a49d70da2572448605250", dependencies = FALSE)
  Install("PredictiveEcology/LandR@a2df2a33fcde78ee16f73565102716fa58917d8b", dependencies = FALSE)
  Install("CeresBarros/ToolsCB@8cdcc5494fdb48c3a3df47c93b2d2cc65c21ce96", dependencies = FALSE)
}


## RESTART R AND KEEP GOING
pkgDir <- file.path(
  if (Sys.info()[["user"]] %in% c("rstudio", "root")) "packages_docker" else "packages",
  version$platform,
  paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
)

.libPaths(pkgDir)
