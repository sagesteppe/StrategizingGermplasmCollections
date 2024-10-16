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
      'Argument to `projection` is required. A good choice for North America is XXXX.')
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

#' Create a regular tessellating grid over the geographic range
#' 
#' This function creates `n` grid cells over a geographic area (`x`), typically a species. It is intended to deal with 'holes' in species ranges to deliver an adequate number of `n`. 
#' @param x An SF object or terra spatraster. the range over which to generate the clusters.
#' @param n Numeric. the number of grid cells desired. Defaults to 20. 
GridSample <- function(x){
  
  
}




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







### adding grid

library(tidyverse)
library(sf)
library(spData)

#' It appears that the default values from sf are adequate to minimize the number of grids, and to decrease the variance in their size over the target area. 
#' 
#' @param x a species range. 
testGridSizes <- function(x){
  
  bound <- st_bbox(target)
  x_dist <- bound['xmax'] - bound['xmin']
  y_dist <- bound['ymax'] - bound['ymin']
  
  ratio <- x_dist/y_dist
  rm(bound, x_dist, y_dist)
  # values < 0.8 indicate y is considerable greater than x
  # values near 1 indicate a natural symmetry between x and y, both values can start at 5. 
  # values < 0.9 > 0.8 indicate y is greater than x
  # values > 1.4 indicate x is much greater than y
  
  if(ratio < 0.8){ # these areas are very long
    x_start = 4; y_start = 7} else if(ratio > 0.8 & ratio < 0.9) {
      x_start = 4; y_start = 6} else if(ratio > 0.9 & ratio < 1.2) { # equilateral
        x_start = 5; y_start = 5} else {
          x_start = 6; y_start = 4} 
  
  gr <- sf::st_make_grid(target, n = c(x_start, y_start), square = FALSE)
  gr <- sf::st_intersection(gr, target) 
  areas <- as.numeric(sf::st_area(gr))
  areas <- sort(areas / max(areas) * 100, decreasing = TRUE)[1:20]
  var_Original <- var(areas,  na.rm = TRUE)
  
  # try with bigger grids
  gr_larger <- sf::st_make_grid(target, n = c(x_start-1, y_start-1), square = FALSE)
  gr_larger <- sf::st_intersection(gr_larger, target) 
  gr_larger_area <- as.numeric(sf::st_area(gr_larger))
  gr_larger_area <- sort(gr_larger_area / max(gr_larger_area) * 100, decreasing = TRUE)[1:20]
  var_larger <- var(gr_larger_area, na.rm = TRUE)
  
  # try with biggest grids
  gr_largest <- sf::st_make_grid(target, n = c(x_start-2, y_start-2), square = FALSE)
  gr_largest <- sf::st_intersection(gr_largest, target) 
  gr_largest_area <- as.numeric(sf::st_area(gr_largest))
  gr_largest_area <- sort(gr_largest_area / max(gr_largest_area) * 100, decreasing = TRUE)[1:20]
  var_largest <- var(gr_largest_area, na.rm = TRUE)
  
  # try with smaller grids
  gr_smaller <- sf::st_make_grid(target, n = c(x_start+1, y_start+1), square = FALSE)
  gr_smaller <- sf::st_intersection(gr_smaller, target) 
  gr_smaller_area <- as.numeric(sf::st_area(gr_smaller))
  gr_smaller_area <- sort(gr_smaller_area / max(gr_smaller_area) * 100, decreasing = TRUE)[1:20]
  var_smaller <- var(gr_smaller_area, na.rm = TRUE)
  
  # try with smallest grids
  gr_smallest <- sf::st_make_grid(target, n = c(x_start+2, y_start+2), square = FALSE)
  gr_smallest <- sf::st_intersection(gr_smallest, target) 
  gr_smallest_area <- as.numeric(sf::st_area(gr_smallest))
  gr_smallest_area <- sort(gr_smallest_area / max(gr_smallest_area) * 100, decreasing = TRUE)[1:20]
  var_smallest <- var(gr_smallest_area, na.rm = TRUE)
  
  results <- data.frame(
    Name = c('Smallest', 'Smaller', 'Original',   'Larger', 'Largest'),
    Grids = c(
      length(gr_smallest), length(gr_smaller), 
      length(gr), length(gr_larger), length(gr_largest)), 
    Variance = c(var_smallest, var_smaller, var_Original, var_larger, var_largest),
    GridNOx = c(x_start+2, x_start+1, x_start, x_start-1, x_start-2), 
    GridNOy = c(y_start+2, y_start+1, y_start, y_start-1, y_start-2)
  )
  return(results)
}

target <- spData::us_states |> 
  dplyr::filter(NAME == 'Rhode Island') |>
  sf::st_transform(32615)

out <- testGridSizes(target)
plot(out$Grids, out$Variance)
abline(v=20, col="blue")

# select the grid size within the acceptable range of grid cells which 
# decrease variance the most. 

# - how do chooose the acceptable number of grids... 

library(sf)
library(tidyverse)
# Determine the size of each grid. 
gr <- sf::st_make_grid(target, n = c(5, 5), square = FALSE) # toy data. 
gr <- sf::st_intersection(gr, target) 
gr <- sf::st_collection_extract(gr, 'POLYGON')

grid_areas <- sf::st_as_sf(gr) |> 
  dplyr::mutate(
    ID   = 1:dplyr::n(),
    Area = as.numeric(sf::st_area(gr))
    )

# order polygons by size, all polygons > 20 will be merged with a neighboring polygon
indices <- grid_areas$Area >= sort(grid_areas$Area, decreasing = TRUE)[20]
to_merge_sf <- gr[!indices,]
merge_into_sf <- gr[indices,] |> sf::st_as_sf()
# before calculating neighbors, we will union, adjacent polygons which we will 
# end up merging to the polygons we will keep. This SHOULD allow for better 
# distribution of there areas into the remaining polygons, making the kept polygons more equal in size. 

gr <- to_merge_sf |> 
  sf::st_union() |> 
  sf::st_cast('POLYGON') |> 
  sf::st_as_sf()  %>% # gotta use pipe to set position of tibbles 
  dplyr::bind_rows(
    gr[indices,] |> sf::st_as_sf(), .
  )

# Determine neighboring polygons
neighbors <- spdep::poly2nb(gr, queen = FALSE)[21:nrow(gr)]

grid_areas <- sf::st_as_sf(gr) |> 
  mutate(
    ID   = 1:dplyr::n(),
    Area = as.numeric(sf::st_area(gr))
  )

full_sized_neighbors <- which( # consider these to be full sized grids
  grid_areas$Area[1:20] / max(grid_areas$Area) >= 0.975) 

# identify neighboring polygons
to_merge_sf <- gr[21:nrow(gr),]
merge_into_sf <- gr[1:20,] 

ggplot() + 
  geom_sf(data = to_merge_sf, fill = 'green') + 
  geom_sf(data = merge_into_sf, fill = 'red')


area2be_reassigned <- vector(length = length(neighbors))
areas <- vector(mode = 'list', length = length(neighbors))
prop_areas <- vector(mode = 'list', length = length(neighbors))
prop_donor <- numeric(length = length(neighbors))

for (i in seq_along(neighbors)){
  
  area2be_reassigned[i] <- sf::st_area(to_merge_sf[i,])
  areas[[i]] <- as.numeric(sf::st_area(gr[neighbors[[i]],]))
  
  prop_donor[i] <- area2be_reassigned[i] / (sum(areas[[i]]) + area2be_reassigned[i]) 
  prop_areas[[i]] <- areas[[i]] / (sum(areas[[i]]) + area2be_reassigned[i]) 
}

rm(area2be_reassigned)

prop_target <- vector(mode = 'list', length = length(prop_areas))
area_sort <- vector(mode = 'list', length = length(prop_areas))
area_des <- vector(mode = 'list', length = length(prop_areas))
nf_pct <- vector(mode = 'list', length = length(prop_areas))
props <- vector(mode = 'list', length = length(prop_areas))
# Using the polygons which will be merged, try to make the following polygons
# as equally sized as possible - without ever removing area from an existing grid. 
for (i in seq_along(area_sort)){
  
  area_des <- (sum(prop_areas[[i]]) + prop_donor[i]) / length(prop_areas[[i]])
  
  if(all(prop_areas[[i]] < area_des)==TRUE){

    # these polygons will all be the same size!... roughly... 
    prop_target[[i]] <- rep(area_des, times = length(prop_areas[[i]]))
    
  } else if(any(prop_areas[[i]] > area_des)==TRUE){
    
    prop_target[[i]] <- numeric(length(prop_areas[[i]]))
    kp <- prop_areas[[i]] < area_des # determine which grids are to big
    area_des <- (sum(prop_areas[[i]][kp]) + prop_donor[i]) / 
      length(prop_areas[[i]][kp])
    
    # make grids smaller than the goal threshold size the threshold, 
    # return grids larger than the threshold size as they are. 
    prop_target[[i]][kp] <- area_des
    prop_target[[i]][-kp] <-prop_areas[[i]][-kp]
  }
  
  nf_pct[[i]] <- setNames( # the existing cover for each grid. 
    prop_areas[[i]] * 100, 
    neighbors[[i]]
  )
  
  props[[i]] <- setNames( # the desired cover for each grid 
    prop_target[[i]] * 100, 
    neighbors[[i]]
  )
}

rm(area_des, area_sort, i, prop_donor, prop_target, areas)

#' place random points in the polygon which will be dissolved with the larger polygons
assignGrid_pts <- function(neighb_grid, focal_grid, props, nf_pct){
  
  # Ensure that we are only working on a grid with 100 points
  pts <- sf::st_sample(focal_grid, size = 100, type = 'regular') |>
    sf::st_as_sf() |> 
    dplyr::mutate(ID = 1:dplyr::n())
  
  samp <- 100
  while (nrow(pts) < 100) {
    if (nrow(pts) < 100) {
      pts <- sf::st_sample(focal_grid, size = samp, type = 'regular') |>
        sf::st_as_sf() |> 
        mutate(ID = 1:n())
    }
    samp <- samp + 1
  } 
  pts <- pts[sample(1:nrow(pts), size = 100, replace = F), ]
  rm(samp)
  
  # identify the nearest neighbor which the points can be assigned to. 
  # we use these to determine what are the 'neediest' neighbors and assign them 
  # points first 
  pts$nf <- sf::st_nearest_feature(pts, neighb_grid) 
  pts$nf <- sf::st_drop_geometry(neighb_grid)[pts$nf,'NAME']
  
  # each grid now receives either the maximum number of nearest neighbors if
  # the desired proportion is lower than the existing nearest neighbors, 
  # or the proportion of nearest points meeting the desired proportion
  
  dists <- sf::st_distance(pts, neighb_grid)
  dists <- data.frame(apply(dists, 2, as.numeric))
  colnames(dists) <- neighb_grid$NAME
  
  assign_pts_frst <- function(x, props, nf_pct){
    
    # ensure these are in the same order. so we can match by position
    # in the loop. 
    need_most2least <- names(sort(nf_pct - props))
    props <- props[need_most2least]
    x <- x[,need_most2least]
    
    x$Assignment  <- NA; x$ID <- 1:nrow(x)
    for(i in 1:length(props)){
      
      # we assign the grids to the neediest grids first, and then work back 
      # removing these points so they are not overwritten. 
      x_sub <- x[is.na(x$Assignment),]
      indices <- x_sub[sort(x_sub[,i], index.return = TRUE)$ix [1:props[i]],'ID']
      x[indices, 'Assignment'] <- names(props)[i]
    }
    return(x)
  }
  
  # if only one grid remains, assign all remaining points to it. 
  frst_assignments <- assign_pts_frst(dists, props = props, nf_pct)
  pts$Assigned <- frst_assignments$Assignment
  rm(nf_pct, dists)
  
  if(any(is.na(frst_assignments$Assignment))){
    
    needAssigned <- pts[is.na(pts$Assigned),]
    dist_final_pts <- data.frame(
      matrix(
        t(sf::st_distance(needAssigned, pts)), ncol = nrow(needAssigned)
      )
    )
    dist_final_pts[dist_final_pts==0] <- NA
    dist_final_pts$ID <- pts$ID
    
    ob <- vector(mode = 'list', length = ncol(dist_final_pts)-1)
    for (i in 1:ncol(dist_final_pts)){
      ob[[i]] <- dist_final_pts[order(dist_final_pts[,i]),  'ID']
    }
  }
  
  rm(frst_assignments)
  # a few points may remain in the center of the object. 
  # We will now try to assign all of these to the groups which do
  # have not yet come close to adequate representation. 
  
  # we will determine what the closest neighbors to the unassigned points are. 
  # to grab them we will buffer the unassigned point and pull out the close neighbors. 
  if(exists('needAssigned')){
    nn <- spdep::knearneigh(pts, k=4)[['nn']][needAssigned$ID,]
    
    for (i in 1:nrow(needAssigned)){
      needAssigned[i, 'Assigned'] <- 
        names(
          which.max(
            table(
              sf::st_drop_geometry(pts)[
                unlist(
                  if(length(nn)==4){nn[1:4]} else {nn[i,1:4]}),
                'Assigned'])
          )
        ) 
    }
    pts <- dplyr::filter(pts, ! ID %in% needAssigned$ID) |>
      dplyr::bind_rows(needAssigned) |>
      dplyr::select(Assigned, ID, geometry = x)
  }
  
  if(exists('needAssigned')){rm(frst_assignments)}
  # Determine if there are points which are 'disconnected' from their remaining neighbors
  nn <- spdep::knearneigh(pts, k=4)[['nn']]
  
  indices <- vector(mode = 'list', length = nrow(pts))
  neighs <- vector(mode = 'list', length = nrow(pts))
  focal <- vector(mode = 'list', length = nrow(pts))
  matches <- vector( length = nrow(pts))
  for (i in 1:nrow(pts)){
    
    indices[[i]] <- nn[i,]
    neighs[[i]] <- 
      names(table(sf::st_drop_geometry(pts)[indices[[i]],'Assigned']))
    focal[[i]] <- sf::st_drop_geometry(pts)[i,'Assigned']
    matches[i] <- focal[[i]] %in% neighs[[i]]
    
  }
  
  rm(focal, indices, nn)
  # if a point is disconnected from it's remaining neighbors then assign it to the neighbor which needs more points to reach the ideal sample ratio for it's grid class.
  if(any(matches==F)){
    indices <- neighs[which(matches==F)] 
    
    realized <- table(pts$Assigned)/ nrow(pts) *100
    props [!names(props) %in% names(realized)] <- 0 # if a small grid is missing say so
    diff <- realized - props # the first entry below will gain the point. 
    
    for (i in 1:length(indices)){
      pts[i,'Assigned'] <-
        names(sort(diff[names(diff) %in% unlist(neighs[i])], decreasing = FALSE)[1])
    }
  }
  
  
  pts <- pts |>
    dplyr::select(Assigned, geometry = x)
  return(pts)
}

neighb_grid <- spData::us_states |> 
  dplyr::filter(
    NAME %in% 
      c('California', 'Oregon', 'Idaho', 'Utah',  'Arizona')) |>
  sf::st_transform(32617)

focal_grid <- spData::us_states |>
  filter(NAME == 'Nevada')|>
  sf::st_transform(32617)

to_merge_sf <- split(to_merge_sf, f = 1:nrow(to_merge_sf))
test <- assignGrid_pts(
  neighb_grid = neighbors, 
  focal_grid = to_merge_sf, 
  props = props, 
  nf_pct = nf_pct)

#' turn the point grid suggestions made by assignGrid_pts into polygons
#' 
#' @param x output of `assignGrid_pts`
#' @param neighb_grid the neighboring grid options
#' @param focal_grid the grid to reassign the area of. 

snapGrids <- function(x, neighb_grid, focal_grid){
  
  # combine the output of `assignGrid_pts` together and create a small
  # polygon around there extent. 
  x <- x |>
    dplyr::group_by(Assigned) |> 
    summarise(geometry = sf::st_combine(geometry)) |> 
    sf::st_convex_hull() |> 
    sf::st_difference() |> 
    sf::st_intersection(focal_grid)
  
  # now place points all across the grid to be divided. 
  pts <- sf::st_sample(focal_grid, 10000) |>
    sf::st_as_sf()
  
  #assign each of the newly generated points to the nearest polygon 
  # created above. The closest polygon will become the points
  # identity. 
  pts$Assigned <- sf::st_drop_geometry(
    x$Assigned)[ sf::st_nearest_feature(pts, x)] 
  pts <- pts |>
    dplyr::group_by(Assigned) |>
    dplyr::summarise(geometry = sf::st_union(x)) |>
    sf::st_convex_hull()
  
  # these points were cast to polygons, and now we remove overlapping areas
  all_pts_surface <- sf::st_difference(pts)
  
  # we will determine snap and interval distances by drawing points which should land into the remaining slivers. This will remove any major
  # holes left by the above process. 
  slivers <- rmapshaper::ms_erase(focal_grid, all_pts_surface)
  sliver_pts <- sf::st_sample(slivers, size = 250)
  snap_dist <- ceiling(
    as.numeric(
      max(
        sf::st_distance(
          sliver_pts, 
          all_pts_surface[
            sf::st_nearest_feature(sliver_pts, all_pts_surface),], 
          by_element = TRUE)
      )
    )
  )
  
  # now fill in the gaps across the spatial data set. 
  all_pts_surface <- rmapshaper::ms_simplify(
    all_pts_surface,  keep_shapes = TRUE, snap = TRUE,
    keep = 1, weighting = 0, snap_interval = snap_dist) |>
    sf::st_difference() |>
    sf::st_buffer(snap_dist) |>
    sf::st_difference()
  
  # join the target grids onto this output. 
  final_grids <- neighb_grid |>
    dplyr::select(Assigned = NAME) |>
    dplyr::bind_rows(int) |>
    dplyr::group_by(Assigned) |>
    dplyr::summarize(geometry = sf::st_union(geometry))
  
  # remove small internal holes which may arise from when the
  # created geometries were joined back to the original grid cells. 
  final_grids <- nngeo::st_remove_holes(final_grids, max_area = 1000)
}

out <- snapGrids(x = test, neighb_grid, focal_grid)

ggplot() + 
  geom_sf(data = out)
