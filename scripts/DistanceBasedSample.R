library(gdistance)

Europe <- raster(system.file("external/Europe.grd", package = "gdistance"))

data(genDist)
data(popCoord)

pC <- as.matrix(popCoord[c("x", "y")])

geoDist <- pointDistance(pC, longlat = TRUE)
geoDist <- as.dist(geoDist)
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
plot(dis_clusters)
clusterCut <- cutree(clusters, 7)

clusters <- hclust(w_dist, method = 'ward.D2')
clusterCut <- cutree(clusters, n)

res_clusters <- hclust(resDist, method = 'ward.D2')
plot(res_clusters)

