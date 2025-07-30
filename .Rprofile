## set CRAN repos; use binary linux packages if on Ubuntu
local({
  options("repos" = c(CRAN = "http://cran.rstudio.com",
                      PE = "http://predictiveecology.r-universe.dev/")
          )

  if (Sys.info()[["sysname"]] == "Linux" && grepl("Ubuntu", utils::osVersion)) {
    .os.version <- strsplit(system("lsb_release -c", intern = TRUE), ":\t")[[1]][[2]]
    # options(repos = c(CRAN = paste0("https://packagemanager.rstudio.com/all/__linux__/",
    #                                 .os.version, "/latest")))
  }
})

## package installation location
rver <- paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
pkgDir <- file.path(
  if (Sys.info()[["sysname"]] == "Linux" && rver == "4.1") "packages_docker" else "packages",
  version$platform,
  paste0(version$major, ".", strsplit(version$minor, "[.]")[[1]][1])
)


if (!dir.exists(pkgDir)) {
  dir.create(pkgDir, recursive = TRUE)
}
.libPaths(pkgDir)

