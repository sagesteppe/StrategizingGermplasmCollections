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



nc <- sf::st_read(system.file("shape/nc.shp", package="sf")) |>
  dplyr::select(NAME) |>
  st_transform(32617)
set.seed(41)
existing_collections <- nc[sample(5),] |>
  sf::st_point_on_surface()

ggplot() + 
  geom_sf(data = existing_collections, color = 'red') + 
  geom_sf(data = zones, alpha = 0.5) + 
  geom_sf(data = kmeans_centers)

zones <- st_voronoi(st_union(existing_collections), st_union(nc)) |>
  st_collection_extract(type = "POLYGON") |> 
  st_sf() |> 
  st_intersection(st_union(nc)) 
zones <- st_intersection(zones, st_union(nc))


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
florida <- spData::us_states |> 
  dplyr::filter(NAME == 'Florida') |>
  sf::st_transform(32617)

bound <- st_bbox(florida)
x_dist <- bound['xmax'] - bound['xmin']
y_dist <- bound['ymax'] - bound['ymin']
if(x_dist > y_dist){x_start = 5; y_start = 4} else {y_start = 5; x_start = 4}

gr <- st_make_grid(florida, n = c(x_start, y_start), square = FALSE)
ints <- sum(lengths( st_intersects(gr, florida) ) > 0) # 12  grids s

gr6 <- st_make_grid(florida, n = c(6, 6), square = FALSE)
gr6 <- gr6[ lengths( st_overlaps(gr6, florida) ) > 0, ]
gr6 <- st_intersection(gr6, florida)

ggplot() + 
  geom_sf(data = florida) + 
  # geom_sf(data = gr, fill = NA) + 
  # geom_sf(data = gr2, fill = NA) + 
  geom_sf(data = gr6, fill = NA) 

# Determine the size of each grid. 

gr6 <- st_as_sf(gr6) |> 
  mutate(Area = as.numeric(st_area(gr6)))

polys2merge <- arrange(gr6, Area) |>
  head(n = nrow(gr6) - 20) 

# Determine neighboring polygons

neighbors <- spdep::poly2nb(gr6, queen = FALSE)
coords <- st_coordinates(st_centroid(st_geometry(gr6)))
plot(st_geometry(gr6), border = "grey")
plot(neighbors, coords, add = TRUE)

# order polysgons by size, all polygons > 20 will be merged with a neighboring polygon
indices <- gr6$Area >= sort(gr6$Area, decreasing = TRUE)[20] 
gr6_neigh2snap2 <- neighbors[indices]
gr6_merge_into <- gr6[indices, ]

full_sized_neighbors <- which(gr6$Area / max(gr6$Area) >= 0.975) # consider these to be full sized grids

# identify neighboring polygons
gr6_neighs2remove <- neighbors[!indices]

# identify the relative sizes of the neighboring polygons

size_props <- function(x, data, full_sized_neighbors){
  
  # if there are neighbors which are full grid size, they will get no records IF
  # there is at least one or more smaller neighbor grid. 
  if(length(x) > 1){x_sub <- x[!x %in% full_sized_neighbors]}
  if(exists('x_sub') == TRUE){
    if(length(x_sub > 1)){x <- x_sub} else {x <- x}
  }
  
  grids <- sf::st_drop_geometry(data)
  areas <- grids[x, 'Area']
  totalArea <- sum(as.numeric(areas))
  
  # if multiple polygon neighbors exist, determine their sizes relative to each other
  recs2receive <- totalArea / (as.numeric(areas)) * 10 # percent of records to receive
  replace(recs2receive, recs2receive == 10, 100)
}

o <- lapply(
  gr6_neighs2remove, 
  FUN = size_props, 
  data = gr6, 
  full_sized_neighbors = full_sized_neighbors)

o

# place random points in the polygon which will be dissolved with the larger polygons





florida <- spData::us_states |> 
  dplyr::filter(
    NAME %in% 
      c('Kansas', 'Wyoming', 'Nebraska', 'Utah',  'New Mexico', 'Oklahoma')) |>
  sf::st_transform(32617)

co <- spData::us_states |>
  filter(NAME == 'Colorado')|>
  sf::st_transform(32617)

pts <- st_sample(co, size = 100, type = 'regular') |>
  st_as_sf() |> 
  mutate(ID = 1:n())

# identify the nearest neighbor which the points can be assigned to. 
pts$nf <- sf::st_nearest_feature(pts, florida) 
pts$nf <- sf::st_drop_geometry(florida)[pts$nf,'NAME']

nrow(pts) / table(pts$nf) # percent nearest feature. 
nf_pct <- setNames(
  as.numeric(table(pts$nf)) / nrow(pts) * 100, 
  names(table(pts$nf))
)

props <- setNames(
  c(15, 10, 30, 2.5, 17.5, 25), 
  names(nf_pct)
)

# each grid now receives either the maximum number of nearest neighbors if
# the desired proportion is lower than the existing nearest neighbors, 
# or the proportion of nearest points meeting the desired proportion

d <- st_distance(pts, florida)
d <- data.frame(apply(d, 2, as.numeric))
colnames(d) <- florida$NAME


assign_pts_frst <- function(x, props, nf_pct){
  
  # ensure these are in the same order. so we can match by position
  # in the for loop. 
  need_most2least <- names(sort(nf_pct - props))
  props <- props[need_most2least] ; x <- x[,need_most2least]
  
  x$Assignment  <- NA
  x$ID <- 1:nrow(x)
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
out <- assign_pts_frst(d, props = props, nf_pct)
pts$ASS <- out$Assignment

if(any(is.na(out$Assignment))){
  
  needAssigned <- pts[is.na(pts$ASS),]
  dist_final_pts <- data.frame(
    matrix(
      t(sf::st_distance(needAssigned, pts)), ncol = 2
    )
  )
  dist_final_pts[dist_final_pts==0] <- NA
  dist_final_pts$ID <- pts$ID
  
  ob <- vector(mode = 'list', length = ncol(dist_final_pts)-1)
  for (i in 1:ncol(dist_final_pts)){
     ob[[i]] <- dist_final_pts[order(dist_final_pts[,i]),  'ID']
  }
  
}

# points are likely to remain in the center of the object. 
# We will now try to assign all of these to the groups which do
# have not yet come close to adequate representation. 

# we will determine what the closest neighbors to the unassigned points are. 
# to grab them we will buffer the unassigned point and pull out the close neighbors. 
dat_mat <- sf::st_distance(pts) |> 
  as.data.frame() |>
  apply(2, as.numeric)
dat_mat[dat_mat==0] <- NA
buf_dist <- ceiling(mean(apply(dat_mat, 2, min, na.rm=TRUE)) * 0.001) * 1000
buffered <- sf::st_buffer(needAssigned, dist = buf_dist)
intersects <- sf::st_intersects(buffered, pts)

nn <- spdep::knearneigh(pts, k=4)[['nn']][needAssigned$ID,]

for (i in 1:nrow(needAssigned)){
  needAssigned[i, 'ASS'] <- 
  names(
    which.max(
      table(sf::st_drop_geometry(pts)[unlist(intersects[[i]]),'ASS'], exclude = NA)
      )
    )
}

pts <- dplyr::filter(pts, ! ID %in% needAssigned$ID) |>
  dplyr::bind_rows( needAssigned ) |>
  dplyr::select(Assigned = ASS, geometry = x)


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

# if a point is disconnected from it's remaining neighbors then assign it to the neighbor which needs more points to reach the ideal sample ratio for it's grid class.  
if(any(matches==F)){
  indices <- neighs[which(matches==F)] 
  
  realized <- table(pts$Assigned)/ nrow(pts) *100
  diff <- realized - props # the first entry below will gain the point. 
  
  for (i in 1:length(indices)){
    pts[i,'Assigned'] <-
      names(sort(diff[names(diff) %in% unlist(neighs[i])], decreasing = FALSE)[1])
  }
}



