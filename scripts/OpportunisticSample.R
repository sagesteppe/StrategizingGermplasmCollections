#' Design additional collections around already existing collections

library(sf)
library(tidyverse)
nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
  dplyr::select(NAME) |>
  sf::st_transform(32617)
set.seed(41)
existing_collections <- nc[sample(5),] |>
  sf::st_point_on_surface()

ggplot() + 
  geom_sf(data = existing_collections, color = 'red') + 
  geom_sf(data = zones, alpha = 0.5) + 
  geom_sf(data = kmeans_centers)

zones <- sf::st_voronoi(st_union(existing_collections), st_union(nc)) |>
  sf::st_collection_extract(type = "POLYGON") |> 
  sf::st_sf() |> 
  sf::st_intersection(st_union(nc)) 
zones <- sf::st_intersection(zones, sf::st_union(nc))


reg_pts <- sf::st_sample(nc, size = 500, type = 'regular', by_polygon = F) |>
  sf::st_coordinates(pts)
fixed_pts <- matrix(rep(sf::st_coordinates(existing_collections), each = 100), ncol = 2)
pts <- rbind(reg_pts, fixed_pts)

kmeans_res <- kmeans(pts, centers = 20)
pts$Cluster <- kmeans_res$cluster

kmeans_centers <- setNames(
  data.frame(kmeans_res['centers'], 1:nrow(kmeans_res['centers'][[1]])), 
  # use the centers as voronoi cores ... ?
  c('X', 'Y', 'Cluster'))
#kmeans_centers_10 <- sf::st_as_sf(kmeans_centers, coords = c('X', 'Y'), crs = 32617)
kmeans_centers_100 <- sf::st_as_sf(kmeans_centers, coords = c('X', 'Y'), crs = 32617)

ggplot() + 
  geom_sf(data = existing_collections, color = 'red') + 
  geom_sf(data = zones, alpha = 0.5) + 
  geom_sf(data = kmeans_centers_10, alpha = 0.5, color = 'blue') + 
  geom_sf(data = kmeans_centers_100, alpha = 0.5, color = 'purple')


