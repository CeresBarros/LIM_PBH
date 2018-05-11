## ---------------------------------
## DAVID ANDISON FIRE DATA
##
## Selection of fires for which further vegetation data can be found
## ---------------------------------

rm(list = ls()); gc(reset = TRUE)

## using - as of April 19th, 2018
# loading reproducible     0.1.4.9015
# loading quickPlot        0.1.3.9002
# loading SpaDES.core      0.1.1.9009
# loading SpaDES.tools     0.1.1.9005
# loading SpaDES.addins    0.1.1
# devtools::install_github("PredictiveEcology/reproducible@development")
# devtools::install_github("PredictiveEcology/quickPlot@development")
# devtools::install_github("PredictiveEcology/SpaDES.core@development")
# devtools::install_github("PredictiveEcology/SpaDES.tools@development")

## requires
library(SpaDES); library(sf); library(ggplot2); library(data.table); library(dplyr)
source("R/R_tools/Useful_functions.R")

## define paths
setPaths(cachePath = file.path("R/SpaDES/cache"),
         modulePath = file.path("R/SpaDES/m"),
         inputPath = file.path("R/SpaDES/inputs"),
         outputPath = file.path("R/SpaDES/outputs"))


## LOAD DATA ---------------------------------------

files = c("albertafires1_postfire", "albertafires2_postfire", "saskatchewanfires_postfire")
folder = "~/../OneDrive/Documents/LandscapesInMotion/data/fires_Dave/Projected_renamed"

# DA_fires <- prepInputs(targetFile = file.path(dataDIR, "albertafires1_postfire.shp"),
#                        url = "https://drive.google.com/file/d/1wGCqI_X1t-PDM5eO6JWQW9hpNv_8zlum/view?usp=sharing",
#                        archive = "Projected_renamed.zip",
#                        alsoExtract = file.path(dataDIR, c("albertafires1_postfire.", "albertafires2_postfire.", 
#                                                           "albertafires1_pretfire.", "albertafires2_prefire.",
#                                                           "saskatchewanfires_postfire.","saskatchewanfires_prefire.")),
#                        fun = "st_read", pkg = "sf",
#                        destinationPath = dataDIR)

for(x in files) {
  eval(parse(text = paste0(
    x, " <- st_read(file.path(folder", ", paste0('", x,"', '.shp')))"
  )))
}

head(albertafires1_postfire)
head(albertafires2_postfire)
head(saskatchewanfires_postfire)

## Saskatchewan attributes table has some repeated columns, check if they are duplicated before deleting
if(all(saskatchewanfires_postfire$FIRE_NAME_ == saskatchewanfires_postfire$FIRE_NAME)) {
  saskatchewanfires_postfire$FIRE_NAME_ <- NULL
}
if(all(saskatchewanfires_postfire$AREA_ == saskatchewanfires_postfire$AREA)) {
  saskatchewanfires_postfire$AREA_ <- NULL
}
if(all(saskatchewanfires_postfire$PERIMETER_ == saskatchewanfires_postfire$PERIMETER)) {
  saskatchewanfires_postfire$PERIMETER_ <- NULL
}

## Uniformize column names
names(albertafires2_postfire)[grep("YEAR", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("EAR", names(saskatchewanfires_postfire))] = 
  grep("YEAR", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("NAME", names(albertafires2_postfire))] = 
  names(saskatchewanfires_postfire)[grep("NAME", names(saskatchewanfires_postfire))] = 
  grep("NAME", names(albertafires1_postfire), value = TRUE)

names(albertafires2_postfire)[grep("FIRE_CODE", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("FIRE_NUM", names(albertafires1_postfire))] = 
  grep("FIRE_ID", names(saskatchewanfires_postfire), value = TRUE)

names(albertafires2_postfire)[grep("ELEMENT", names(albertafires2_postfire))] = 
  names(albertafires1_postfire)[grep("SURVIVAL", names(albertafires1_postfire))] = 
  names(saskatchewanfires_postfire)[grep("CLASS", names(saskatchewanfires_postfire))] = "SEV_CLAS"


ABSK_firePerims <- rbind(albertafires1_postfire[, "FIRE_NAME"], albertafires2_postfire[, "FIRE_NAME"], saskatchewanfires_postfire[, "FIRE_NAME"])
ABSK_firePerims <- ABSK_firePerims %>% group_by(FIRE_NAME) %>% 
  summarise(Group = unique(FIRE_NAME)) %>% 
  st_cast()

# dev(); plot(ABSK_firePerims)
CA_admin <- Cache(prepInputs, url = "http://www12.statcan.gc.ca/census-recensement/2011/geo/bound-limit/files-fichiers/gpr_000a11a_e.zip",
                  destinationPath = "data/CA_admin",
                  cacheRepo = "analyses/cache")
CA_admin <- CA_admin[CA_admin$PRENAME %in% c("Alberta", "Saskatchewan"),]
CA_admin <- sp::spTransform(CA_admin, CRSobj = st_crs(ABSK_firePerims)$proj4string)

FMAs <- Cache(prepInputs, targetFile = "FMA_Boundary_Updated.shp",
              url = "https://drive.google.com/file/d/1nTFOcrdMf1hIsxd_yNCSTr8RrYNHHwuc/view?usp=sharing", 
              destinationPath = "data/FMA", studyArea = CA_admin, writeCropped = FALSE,
              cacheRepo = "analyses/cache")

CA_admin <- st_as_sfc(CA_admin)   ## plots do not align well when CA_admin is an sf - doesnt need to be one.
FMAs <- st_as_sf(FMAs)

plot(CA_admin); plot(FMAs, add = TRUE); plot(ABSK_firePerims, add = TRUE)


## fire data with FMAs
FMAs <- st_transform(FMAs, crs = st_crs(ABSK_firePerims))
firesFMAs <- Cache(st_intersection, x = FMAs, y = ABSK_firePerims, cacheRepo = "analyses/cache")

plot(CA_admin); plot(FMAs, add = TRUE); plot(firesFMAs, add = TRUE)

## subset to AB only
AB.poly <- CA_admin[1]
AB.poly <- st_transform(AB.poly, crs = st_crs(firesFMAs))

## group fires by FMA
firesFMAs_AB <- Cache(st_intersection, x = firesFMAs, y = AB.poly)

plot(firesFMAs_AB[, "Name"])
dir.create("analyses/firesInFMAs")
st_write(firesFMAs_AB, dsn = "analyses/firesInFMAs/firesDA_FMAs_AB.shp", delete_dsn = TRUE)
