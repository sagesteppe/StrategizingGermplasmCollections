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



?commuteDistance

# NOW ATTEMPT WITH ANOTHER VARIABLE. #

# manually (arbitrarily) setting cluster number works.... 
resDist_scaled <- dist(scale(resDist), method = 'euclidean') # scale variables
geo_clusters <- hclust(resDist_scaled,  method = 'ward.D2') # create a hierarchical 
# clustering tree. 
plot(geo_clusters)

####### THIS WORKS, BUT WE ARE UNABLE TO USE ALL INDICES. 
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


