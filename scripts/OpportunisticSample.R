#' Design additional collections around already existing collections

library(sf)
library(tidyverse)
nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
  dplyr::select(NAME) |>
  sf::st_transform(32617)
set.seed(42)
existing_collections <- nc[sample(1:nrow(nc),size = 5),] |>
  sf::st_point_on_surface()


#' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param n_col Numeric. The total number of desired collections. Defaults to 20.
#' @param collections an sf point geometry data set of where existing collections have been made.
VoronoiSampler <- function(polygon, n_col, collections){
  
  pts <- sf::st_sample(polygon, size = 20, type = 'regular', exact = TRUE) |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = x)
  
  pts <- dplyr::bind_rows(
    existing_collections, 
    pts[-st_nearest_feature(existing_collections, pts),], ) |>
    dplyr::slice_head(n=20)
  
  vorons <- sf::st_voronoi(sf::st_union(pts), sf::st_as_sfc(st_bbox(polygon)))
  vorons <- sf::st_intersection(st_cast(vorons), sf::st_union(polygon))
  variance <- var(sf::st_area(vorons))
  
  return(list(
    Variance = variance, 
    Polygons = vorons
  ))
}

lapply(seq_len(10), FUN = VoronoiSampler, polygon = nc, n_col = 10, collections = existing_collections)







pts <- st_sample(nc, size = 20, type = 'regular', exact = TRUE) |>
  st_as_sf() |>
  rename(geometry = x)

locs <- bind_rows(
  existing_collections, 
  pts[-st_nearest_feature(existing_collections, pts),], ) |>
  slice_head(n=20)

v <- sf::st_voronoi(sf::st_union(locs), sf::st_as_sfc(st_bbox(nc)))
v <- sf::st_intersection(st_cast(v), sf::st_union(nc))
var(sf::st_area(v))

ggplot() + 
  geom_sf(data = v)  +
  geom_sf(data = existing_collections, color = 'red')

#' get n points located reasonably far away from existing collections
#' 
#' This is a component of finishing up an opportunistic sample design trying to maximize the distance between 
#' the already accessioned collections and new collections. 
#' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param collections an sf point geometry data set of where existing collections have been made. 
FurthestPoint <- function(polygon, collections){
  
  pts <- sf::st_sample(polygon, size = 500, type = 'regular') |>
    sf::st_as_sf()
  d_mat <- sf::st_distance(pts, collections)
  pts$meanDistance <- apply(d_mat, 1, mean)
  x <- sf::st_coordinates(pts)[,1];  y <- sf::st_coordinates(pts)[,2]
  pts <- sf::st_drop_geometry(pts);  pts <- cbind(pts, x, y)
  
  r <- terra::rast(polygon, ncol = 50, nrow = 50)
  mg <- gstat::gstat(
    id = "meanDistance", locations = ~x+y, formula = meanDistance~1,  
    data=pts, nmax=5, set=list(idp = .5))
  z <- terra::interpolate(r, mg, debug.level=0, index=1)
  z <- terra::mask(z, polygon)
  
  v <- terra::as.polygons(z) |>
    sf::st_as_sf()
  
  quant975 <- quantile(v$meanDistance.pred, probs = 0.975)
  # identify the points in the most disconnected distances. 
  polys <- dplyr::mutate(v, 
                         Target = if_else(meanDistance.pred > quant975, TRUE, FALSE), 
                         .before = geometry) |>
    dplyr::filter(Target == TRUE) |>
    sf::st_union() |>
    sf::st_cast('POLYGON') 
  
  if(length(polys)>1){
    polys <- polys[order(sf::st_area(polys), decreasing = TRUE)[1],]
  }
  
  pt <- sf::st_point_on_surface(polys) |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = x) 
  existing_collections <- dplyr::bind_rows(collections, pt)
  
  return(existing_collections)
  
}

ggplot() + 
  geom_sf(data = nc) + 
  geom_sf(data = pt) + 
  geom_sf(data = existing_collections, color = 'red')
  

#' @param x the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param collections an sf point geometry data set of where existing collections have been made. 
#' @param desired_collections the number of collections desired
OpporunisticSample <- function(x, collections, desired_collections){
  
  while(nrow(collections) < desired_collections){
    collections <- FurthestPoint(nc, collections = collections)
  }
  
  return(collections)
  
}

out <- OpporunisticSample(nc, collections = existing_collections, desired_collections = 20)

ggplot() +
  geom_sf(data = nc) + 
  geom_sf(data = out, aes(color = NAME))

install.packages('spsurvey')
library(spsurvey)

legacy <- grts(nc, n_base = 20, legacy_sites = existing_collections)

ggplot() + 
  geom_sf(data = nc) + 
  geom_sf(data = legacy[['sites_legacy']], color = 'red') + 
  geom_sf(data = legacy[['sites_base']])
