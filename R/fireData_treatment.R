## ---------------------------------
## DAVID ANDISON FIRE DATA
##
## Data treatment
## ---------------------------------
library(sf)
library(quickPlot)

data.folder <- "data/fires_Dave/Projected_renamed/"

file.ls <- list.files(data.folder, pattern = "shp")
file.ls <- file.ls[!grepl("water", file.ls)]

for(x in file.ls){
  eval(expr = parse(text = paste0(sub(".shp", "", x), " <- st_read('", file.path(data.folder, x), "')"))) 
}

## subset pre-fire shapefiles to post-fire shapefiles
dev()
Plot(albertafires1_prefire[1])
Plot(albertafires1_postfire[1])

pre_postfire_intersect <- st_intersection(albertafires1_prefire, albertafires1_postfire)
Plot(pre_postfire_intersect[1])


## match the rown.names to those in post-fire data
row.names(pre_postfire_intersect) <- gsub(".* ", "", row.names(pre_postfire_intersect))

keep <- row.names(pre_postfire_intersect)

## now retrive the post fire data that matches the rows that intersected with pre-fire data
## because the same fire polygon intersects with different vegetation polygons, the row names
## appear like "decimals" polyID.x, where x is an index and polyID is the repeated polygon 
postfire.dt <- albertafires1_postfire@data[keep, ]

## change rownames to something sensible
rownames(postfire.dt) <- as.character(1:nrow(postfire.dt))

## Make polygon IDs match the new row names
pre_postfire_intersect <- spChFIDs(pre_postfire_intersect, rownames(postfire.dt))

## finally associate the dataframe with the polygons, by creating a SpatialPolygonsDataFrame - didnt work... only post fire data here...
pre_postfire_intersect <- SpatialPolygonsDataFrame(pre_postfire_intersect, postfire.dt)
