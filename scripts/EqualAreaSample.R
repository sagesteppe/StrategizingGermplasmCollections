#' Create equal area polygons over a geographic range
#' 
#' This function creates `n` geographic clusters over a geographic area (`x`), typically a species
#' range, using kmeans clustering. 
#' @param x An SF object or terra spatraster. the range over which to generate the clusters.
#' @param n Numeric. the number of clusters desired. Defaults to 20. 
#' @param pts Numeric. the number of points to use for generating the clusters, these will be placed in a grid like fashion across `x`. The exact number of points used may deviate slightly from the user submitted value to allow for equidistant spacing across `x`. Defaults to 10,000
#' @param projection Numeric. An EPSG code for a planar coordinate projection, in meters, for use with the function. For species with very narrow ranges a UTM zone may be best (e.g. 32611 for WGS84 zone 11 north, or 29611 for NAD83 zone 11 north). Otherwise a continential scale projection like 5463. See LINK for more information on CRS. 
#' @param returnProjected. Boolean. Whether to return the data set in the original input CRS (FALSE), or in the new `projection` (True). Defaults to FALSE
#' @examples \donttest{
#' nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
#' dplyr::select(NAME)
#'
#' zones <- EqualAreaSample(nc, n = 20, pts = 1000, projection = 32617)
#'
#' plot(nc, main = 'Counties of North Carolina')
#' plot(zones, main = 'Clusters')
#' }
#' @export
EqualAreaSample <- function(x, n, pts, projection, returnProjected){
  
  if(missing(n)){n <- 20}; if(missing(pts)){pts <- 10000}
  if(missing(projection)){
    message(
      'Argument to `projection` is required. A suitable choice for all of North America is 5070.')
    }
  if(missing(returnProjected)){returnProjected <- FALSE}
  
  # determine which portions of the STZ are likely to be populated by converting
  # the sdm output to vectors and masking the STZ to this area. 
  pts <- sf::st_sample(x, size = pts, type = 'regular', by_polygon = F)
  
  kmeans_res <- kmeans(sf::st_coordinates(pts), centers = n)
  pts$Cluster <- kmeans_res$cluster
  
  # gather the geographic centers of the polygons. 
  kmeans_centers <- setNames(
    data.frame(kmeans_res['centers'], 1:nrow(kmeans_res['centers'][[1]])), 
    # use the centers as voronoi cores ... ?
    c('X', 'Y', 'Cluster'))
  kmeans_centers <- sf::st_as_sf(kmeans_centers, coords = c('X', 'Y'), crs = 4326)
  
  voronoi_poly <- kmeans_centers |> # create polygons surrounding the 
    sf::st_transform(projection) |> # clusters using voronoi polygons
    sf::st_union() |>  
    sf::st_voronoi() |>  
    sf::st_cast() |>  
    sf::st_as_sf() |>
    sf::st_make_valid() 
  
  lkp <- c(geometry = "x") 
  
  voronoi_poly <- sf::st_intersection(
    # reduce the extent of the voronoi polygons to the area of analysis. 
    voronoi_poly, sf::st_union(
      sf::st_transform(x, 
                       sf::st_crs(voronoi_poly)))
  ) |> 
    # we can assign an arbitrary number. 
    dplyr::mutate(Cluster = 1:nrow(voronoi_poly)) |>
    sf::st_make_valid() |>
    sf::st_as_sf() |>
    dplyr::rename(any_of(lkp))
  
  if(returnProjected == FALSE){
    voronoi_poly<- sf::st_transform(voronoi_poly, sf::st_crs(x))
    # return the object in the original projection
  }
  
  realN <- nrow(pts) # return the true sample size
  return(voronoi_poly)
  
}

library(tidyverse)
nc <- spData::us_states |>
  filter(NAME == 'Rhode Island')

zones <- EqualAreaSample(nc, n = 20, pts = 1000, projection = 32617)

plot(nc, main = 'Counties of North Carolina')
plot(zones, main = 'Clusters')
