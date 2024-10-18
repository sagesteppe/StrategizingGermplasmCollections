#' Design additional collections around already existing collections

library(sf)
library(tidyverse)
nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
  dplyr::select(NAME) |>
  sf::st_transform(32617)
set.seed(42)
existing_collections <- nc[sample(1:nrow(nc),size = 5),] |>
  sf::st_point_on_surface()

min(variance)
npbs[["t0"]] # the statistic of interest 
npbs[["bca"]][["lower"]]
npbs[["bca"]][["upper"]]

#' Generate a sampling grid using points as the input data - including existing collections
#' 
#' This function utilizes a regular, or nearly so in the case of existing collections, grid of points 
#' to develop a sampling scheme or n polygons. 
#' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
#' @param n Numeric. The total number of desired collections. Defaults to 20.
#' @param collections an sf point geometry data set of where existing collections have been made.
PointBasedSample <- function(polygon, n, collections){
  
  #' Recursively grab a named component of a list. 
  #' @param x a list of lists
  #' @param element the quoted name of the list element to extract. 
  get_elements <- function(x, element) { # @ StackOverflow Allan Cameron 
    if(is.list(x))
    {
      if(element %in% names(x)) x[[element]]
      else lapply(x, get_elements, element = element)
    }
  }
  
  #' Make a voronoi sample of an area n times
  #' 
  #' Split an area up into n polygons of roughly equal area, optionally removing 
  #' some of the default points and replacing them with existing collections to 
  #' build the future collections around. 
  #' 
  #' @param polygon the input sf polygon, i.e. species range or administrative unit, where sampling is desired. 
  #' @param n Numeric. The total number of desired collections. Defaults to 20.
  #' @param collections an sf point geometry data set of where existing collections have been made.
  #' @param reps Numeric. The number of times to rerun the voronoi algorithm, the set of polygons with the most similar sizes, as
  #' measured using their variance of areas will be selected. Defaults to 100. 
  VoronoiSampler <- function(polygon, n, collections, reps){
    
    pts <- sf::st_sample(polygon, size = 20, type = 'regular', exact = TRUE) |>
      sf::st_as_sf() |>
      dplyr::rename(geometry = x)
    
    if(exists('collections')){
      pts <- dplyr::bind_rows(
        existing_collections, 
        pts[-sf::st_nearest_feature(existing_collections, pts),], ) |>
        dplyr::slice_head(n = n)
    } else {pts <- dplyr::slice(head(pts, n = n))}
    
    vorons <- sf::st_voronoi(sf::st_union(pts), sf::st_as_sfc(st_bbox(polygon)))
    vorons <- sf::st_intersection(st_cast(vorons), sf::st_union(polygon))
    variance <- var(as.numeric(sf::st_area(vorons))/10000)
    
    # need to define two slots, one for the variance numeric results, and one for the polygons
    
    return(list( # not the way !
      'Variance' = variance,
      'Polygons' = vorons
    ))
  }
  
  
  # we apply the voronoi process a number of replicated times, defaults to 100
  voronoiPolygons <- replicate(
    reps, 
    VoronoiSampler(polygon = nc, n = 10, collections = existing_collections), 
    simplify = FALSE)
  
  get_elements <- function(x, element) { # @ StackOverflow Allan Cameron 
    if(is.list(x))
    {
      if(element %in% names(x)) x[[element]]
      else lapply(x, get_elements, element = element)
    }
  }
  
  # we use variance to determine the configuration of voronoi polygons which have
  # the most equally sized polygons. 
  variance <- unlist(get_elements(voronoiPolygons, 'Variance'))
  selected <- voronoiPolygons[which.min(variance)][[1]]$Polygons
  
  # Determining the 0.1% quantile for the variance in size of the sampling grids. 
  # Using non-parametric approaches, of bootstrap resampling (replicates = 9999) ,
  # with an 95% confidence level. 
  
  # we can show that the polygon arrangement we have chosen is in the top 1000 of
  # options if npbs[["bca"]][["lower"]] > min(variance) == TRUE . 
  # If the above condition is not meet, we can also say that it is less than the estimate
  # npbs[["t0"]] < min(variance)
  npbs <- nptest::np.boot(x = variance, statistic = quantile, 
                          probs = c(0.001), level = 0.95)
  
  
} 
