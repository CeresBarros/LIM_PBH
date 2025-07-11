## set CRAN repos; use binary linux packages if on Ubuntu
local({
  options("repos" = c(CRAN = "https://cran.rstudio.com",
                      PE = "https://predictiveecology.r-universe.dev/")
          )

  if (Sys.info()[["sysname"]] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
    .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
    # options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
    #                                 .os.version, "/latest")))
  }
})

## package installation location
pkgDir <- file.path(
  if (Sys.info()[["user"]] %in% c("rstudio", "root")) "packages_docker" else "packages",
  version$platform,
  paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
)

## settings for for-cast
if (grepl("for-cast", Sys.info()["nodename"])) {
  data.table::setDTthreads(25)
  options(bitmapType="cairo")
}
if (!dir.exists(pkgDir)) {
  dir.create(pkgDir, recursive = TRUE)
}
# .libPaths(pkgDir)
