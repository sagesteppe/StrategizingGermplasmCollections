library(gdistance)

Europe <- raster(system.file("external/Europe.grd", package = "gdistance"))

data(genDist)
data(popCoord)

pC <- as.matrix(popCoord[c("x", "y")])

geoDist <- pointDistance(pC, longlat = TRUE) # Simple geographic distance
# using a great circle between the different locations. 
Europe <- aggregate(Europe, 3) 

clumps = raster::clump(Europe, directions = 8)
clumps[clumps[] != 1] <- NA
Europe = Europe * clumps


tr <- transition(Europe, mean, directions = 8)
trC <- geoCorrection(tr, "c", scl = TRUE)
trR <- geoCorrection(tr, "r", scl = TRUE)

cosDist <- costDistance(trC, pC)
resDist <- commuteDistance(trR, pC)

dis_clusters <- hclust(cosDist, method = 'ward.D2')
popCoord$Cluster <- cutree(dis_clusters, 7)

res_clusters <- hclust(resDist, method = 'ward.D2')
plot(res_clusters)





### Reprex ###

library(gdistance)

data(popCoord)
geoDist <- pointDistance(
  as.matrix(popCoord[c("x", "y")]), 
  longlat = TRUE) # calculate great circle distances between
# locations. 


# manually (arbitrarily) setting cluster number works.... 
geoDist_scaled <- dist(scale(geoDist), method = 'euclidean') # scale variables
geo_clusters <- hclust(geoDist_scaled,  method = 'ward.D2') # create a hierarchical 
# clustering tree. 
plot(geo_clusters)

popCoord$Cluster <- cutree(geo_clusters, 6) # assign clusters to original data. 
plot(popCoord$x, popCoord$y, col = popCoord$Cluster) # visualize results
text(popCoord$x+0.5,  popCoord$y+0.75, labels = popCoord$Population, cex= 0.6)


####### THIS WORKS, BUT WE ARE UNABLE TO USE ALL INDICES. 
NoClusters <- NbClust::NbClust(
  data = as.dist(geoDist),
  diss = geoDist_scaled, distance = NULL,
  min.nc = 2, max.nc = 20, 
  method = 'complete', index = 'silhouette'
)

popCoord$Cluster <- NoClusters$Best.partition # assign clusters to original data. 
plot(popCoord$x, popCoord$y, col = popCoord$Cluster) # visualize results - so many
# clusters impossible for me to tell te colours apart... 
text(popCoord$x+0.5,  popCoord$y+0.75, labels = popCoord$Population, cex= 0.6)
text(popCoord$x-0.5,  popCoord$y-0.75, labels = popCoord$Cluster, cex= 0.5)


# NOW ATTEMPT WITH ANOTHER VARIABLE. #

# manually (arbitrarily) setting cluster number works.... 
resDist_scaled <- dist(scale(resDist), method = 'euclidean') # scale variables
geo_clusters <- hclust(resDist_scaled,  method = 'ward.D2') # create a hierarchical 
# clustering tree. 
plot(geo_clusters)

####### THIS WORKS, BUT WE ARE UNABLE TO USE ALL INDICES, just silhouette plots. 
NoClusters <- NbClust::NbClust(
  data = as.dist(resDist),
  diss = resDist_scaled, distance = NULL,
  min.nc = 2, max.nc = 20, 
  method = 'complete', index = 'silhouette'
)

popCoord$Cluster <- NoClusters$Best.partition # assign clusters to original data. 
plot(popCoord$x, popCoord$y, col = popCoord$Cluster) # visualize results - so many
# clusters impossible for me to tell the colours apart... 
text(popCoord$x+0.5,  popCoord$y+0.75, labels = popCoord$Population, cex= 0.6)
text(popCoord$x-0.5,  popCoord$y-0.75, labels = popCoord$Cluster, cex= 0.5)











#################################################################################
#                         TRANSITION LAYER PLAY                                #

library(terra)
library(gdistance)

template <- terra::rast('../results/SDM/Raster/Bradypus_test.tif')

domain <- terra::crop(
  template$Supplemented, 
  terra::ext(terra::as.polygons(template$Supplemented))
  )

setwd('~/Documents/assoRted/StrategizingGermplasmCollections/scripts')
p2dat <- '~/Documents/Geospatial/'

f <- paste0(
  'GTOPO30arc_Americas/', 
  c('gt30w060n40.tif', 'gt30w100n40.tif', 'gt30w100s10.tif', 'gt30w060s10.tif'))

elev <- terra::mosaic(terra::sprc(paste0(p2dat, f)))
elev <- terra::crop(elev, domain)
slope <- terra::terrain(elev, "slope", neighbors=8)
slope <- terra::resample(slope, domain)
m_val <- terra::global(slope, max, na.rm = TRUE)
slope <- raster::raster(terra::subst(slope, NA, m_val*10))
rm(f)

tr <- transition(slope, mean, directions = 8)
trC <- geoCorrection(tr, "r", scl = TRUE)

pts <- sf::st_sample(
  sf::st_as_sf(terra::as.polygons(template$Supplemented)),
  size = 500) |>
  sf::st_as_sf() |>
  sf::st_coordinates() |>
  as.matrix()


geoDist <- pointDistance(
  pts, 
  longlat = TRUE) # calculate great circle distances between
# locations. 


# manually (arbitrarily) setting cluster number works.... 
geoDist_scaled <- dist(scale(geoDist), method = 'euclidean') # scale variables
geo_clusters <- hclust(geoDist_scaled,  method = 'ward.D2') # create a hierarchical 

plot(geo_clusters)

NoClusters <- NbClust::NbClust(
  data = as.dist(geoDist),
  diss = geoDist_scaled, distance = NULL,
  min.nc = 2, max.nc = 20, 
  method = 'complete', index = 'silhouette'
)


pts <- data.frame(pts)

pts$Cluster <- NoClusters$Best.partition # assign clusters to original data. 
plot(pts$X, pts$Y, col = pts$Cluster) # visualize results - so many
# clusters impossible for me to tell te colours apart... 
text(pts$X-0.5,  pts$Y-0.75, labels = pts$Cluster, cex= 0.5)





test <- commuteDistance(trC, pts)
pts <- data.frame(pts)

# try out on this data set real quick... 

# manually (arbitrarily) setting cluster number works.... 
resDist_scaled <- dist(scale(test), method = 'euclidean') # scale variables
geo_clusters <- hclust(resDist_scaled,  method = 'ward.D2') # create a hierarchical 
plot(geo_clusters)

pts$Cluster <- cutree(geo_clusters, 7) # assign clusters to original data. 
plot(pts$X, pts$Y, col = pts$Cluster) # visualize results
text(pts$X+0.5,  pts$Y+0.75, cex= 0.6)




plot(slope)




# now burn the rivers into a layer as a feature which is highly resistant to movement. 
rivers <- terra::vect(paste0(p2dat, 'major_rivers/MajorRivers.shp'))
rivers <- terra::crop(rivers, slope)
rivers <- sf::st_as_sf(rivers) |>
  sf::st_transform('+proj=cea +lon_0=-423.6767578 +lat_ts=0 +datum=WGS84 +units=m +no_defs') |>
  sf::st_buffer(dist = 1000) |> # this should be the cell size - it will buffer in both directions
  dplyr::mutate(River = 1)  |> # and prevent any diagonal gaps. 
  sf::st_transform( terra::crs(slope)) 

r <- terra::rast(ext(slope), resolution=res(slope))
terra::crs(r) <- terra::crs(slope)
rivers <- terra::rasterize(vect(rivers), r, field = 'River', background = NA)

# repeat the process with lakes. - this basically for the great lakes. 

sf::st_read(paste0(p2dat, 'GLWD_level1/glwd_1.shp'))
lakes <- terra::vect(paste0(p2dat, 'GLWD_level1/glwd_1.shp'))
terra::crs(lakes) <- 'EPSG:4326'
lakes <- terra::crop(lakes, slope)

lakes <- sf::st_as_sf(lakes) |>
  sf::st_transform('+proj=cea +lon_0=-423.6767578 +lat_ts=0 +datum=WGS84 +units=m +no_defs') |>
  dplyr::filter(TYPE == 'Lake') |>
  dplyr::select(geometry) |> 
  dplyr::mutate(Lake = 1)  |> 
  sf::st_transform( terra::crs(slope)) 

r <- terra::rast(ext(slope), resolution=res(slope))
terra::crs(r) <- terra::crs(slope)
lakes <- terra::rasterize(terra::vect(lakes), r, field = 'Lake', background = NA)

impermeable <- terra::app( c(lakes, rivers), mean, na.rm= TRUE)
plot(impermeable)


# terra::writeRaster(impermeable, 'Rivertest.tif', overwrite = TRUE)
rm(rivers, lakes)
