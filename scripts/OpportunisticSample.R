#' Design additional collections around already existing collections

library(sf)
library(tidyverse)
nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
  dplyr::select(NAME) |>
  sf::st_transform(32617)
set.seed(42)
existing_collections <- nc[sample(1:nrow(nc),size = 5),] |>
  sf::st_point_on_surface()

ggplot() + 
  geom_sf(data = nc)  +
  geom_sf(data = existing_collections) 

#' get n points located reasonably far away from existing collections
#' 
#' This is a component of finishing up an opportunistic sample design trying to maximize the distance between 
#' the already accessioned collections and new collections. 
#' @param x the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param collections an sf point geometry data set of where existing collections have been made. 
#' @param i the index of for loop position. 
FurthestPoint <- function(x, collections, i){
  
  pts <- sf::st_sample(x, size = 500, type = 'regular') |>
    sf::st_as_sf()
  d_mat <- sf::st_distance(pts, collections)
  pts$meanDistance <- apply(d_mat, 1, mean)
  x <- sf::st_coordinates(pts)[,1];  y <- sf::st_coordinates(pts)[,2]
  pts <- sf::st_drop_geometry(pts);  pts <- cbind(pts, x, y)
  
  r <- terra::rast(nc, ncol = 50, nrow = 50)
  mg <- gstat::gstat(
    id = "meanDistance", locations = ~x+y, formula = meanDistance~1,  
    data=pts, nmax=5, set=list(idp = .5))
  z <- terra::interpolate(r, mg, debug.level=0, index=1)
  z <- terra::mask(z, nc)
  
  v <- terra::as.polygons(z) |>
    sf::st_as_sf()
  
  quant975 <- quantile(v$meanDistance.pred, probs = 0.975)
  # identify the points in the last 5% most disconnected distance. 
  polys <- dplyr::mutate(v, 
                         Target = if_else(meanDistance.pred > quant975, TRUE, FALSE), 
                         .before = geometry) |>
    dplyr::filter(Target == TRUE) |>
    sf::st_union() |>
    sf::st_cast('POLYGON') 
  
  if(length(polys)>1){
    ord <- order(sf::st_area(polys), decreasing = TRUE)
  }
  
  pt <- sf::st_point_on_surface(polys) |>
    sf::st_as_sf() |>
    dplyr::rename(geometry = x) |>
    dplyr::mutate(SampledOrder = i, .before = geometry) 
  
  return(pt)
  
}

pt <- FurthestPoint(nc, collections = existing_collections, i = 1)


ggplot() + 
  geom_sf(data = nc) + 
  geom_sf(data = pt)
  

