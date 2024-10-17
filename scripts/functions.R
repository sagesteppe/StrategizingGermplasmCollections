library(tidyverse)
library(sf)

#' place random points in the polygon which will be dissolved with the larger polygons
#' 
#' This function is ran within `GridBasedSample` it will place points throughout polygon geometries
#' which should be merged to larger polygons and assign them to neighboring polygons
#' based on how much area we want to grow these polygons too. 
#' @param neighb_grid 
#' @param focal_grid 
#' @param props 
#' @param nf_pct 
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
  
  if(nrow(neighb_grid)==1){
    pts$Assigned <- sf::st_drop_geometry(neighb_grid$ID)} else {
      
      # identify the nearest neighbor which the points can be assigned to. 
      # we use these to determine what are the 'neediest' neighbors and assign them 
      # points first 
      pts$nf <- sf::st_nearest_feature(pts, neighb_grid) 
      pts$nf <- sf::st_drop_geometry(neighb_grid)[pts$nf,'ID']
      
      # each grid now receives either the maximum number of nearest neighbors if
      # the desired proportion is lower than the existing nearest neighbors, 
      # or the proportion of nearest points meeting the desired proportion
      dists <- sf::st_distance(pts, neighb_grid)
      dists <- data.frame(apply(dists, 2, as.numeric))
      colnames(dists) <- neighb_grid$ID
      
      # if only one grid remains, assign all remaining points to it. 
      frst_assignments <- assign_pts_frst(dists, props = props, nf_pct)
      if(exists('frst_assignments')){pts$Assigned <- frst_assignments$Assignment}
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
          dplyr::rename(geometry = x) |>
          dplyr::bind_rows(needAssigned) |>
          dplyr::select(Assigned, ID, geometry = x)
      }
      pts <- pts[! sf::st_is_empty(pts), ]
      
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
    }
  lkup <- c(geometry = "x")
  
  pts <- pts |>
    dplyr::rename(dplyr::any_of(lkup)) |>
    dplyr::select(Assigned, geometry) |>
    dplyr::mutate(Assigned = as.numeric(Assigned))
  return(pts)
}


# NEED TO FIX SEVERAL ITEMS

# 2) Error in spdep::knearneigh(pts, k = 4)[["nn"]][needAssigned$ID, ] : 
# subscript out of bounds

##################################### SAND BOX ###############################




# Problems with spdep::poly2nb(gr, queen = FALSE) with New York
# Error in Error in (function (msg)  : 
# TopologyException: Input geom 1 is invalid: Self-intersection at 1113087.3264981101 4525176.2385835573 with Indiana
target <- spData::us_states |> 
  dplyr::filter(NAME == 'Florida') |>
  sf::st_transform(32615)

TestGridSizes(target)

ggplot() + 
  geom_sf(data = output, aes(fill = Assigned))

output <- GridBasedSample(target)



x <- target
gr <- sf::st_make_grid(x, n = c(6, 6), square = FALSE) # toy data. 
gr <- sf::st_intersection(gr, x) 
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

full_sized_neighbors <- which( # consider these to be full sized grids
  grid_areas$Area[1:20] / max(grid_areas$Area) >= 0.975) 

# identify neighboring polygons
to_merge_sf <- gr[21:nrow(gr),]
merge_into_sf <- gr[1:20,] 

# if their are no neighbors, this implies that the focal grid area is isolated, e.g. an island
# we will use distance based neighbors for that focal area.  
for (i in seq_along(neighbors)){
  if(sum(neighbors[[i]])==0){
    need_distance_neighbs <- which(sum(neighbors[[i]])==0)
    gr_pos <- sf::st_point_on_surface(gr)
    ob <- spdep::knearneigh(gr_pos, k=4)[['nn']]
    neighbors[[i]] <- ob[20+need_distance_neighbs,]
    neighbors[[i]] <- neighbors[[i]][neighbors[[i]]<=20] # remove neighbors which will be dropped. 
  }
}



from <- sf::st_point_on_surface(to_merge_sf[need_distance_neighbs,])
destinations <- merge_into_sf[neighbors[[1]],]

a <- sf::st_union(from, destinations[1,]) |> sf::st_cast('LINESTRING')
b <- sf::st_union(from, destinations[2,]) |> sf::st_cast('LINESTRING')
c <- sf::st_union(from, destinations[3,]) |> sf::st_cast('LINESTRING')

lengths(sf::st_intersects(a, merge_into_sf[destinations,]))
lengths(sf::st_intersects(b, merge_into_sf[destinations,]))
lengths(sf::st_intersects(c, merge_into_sf[destinations,]))

ggplot() + 
  geom_sf(data = gr) + 
  geom_sf(data = a) + geom_sf(data = b) + geom_sf(data = c) +  
  geom_sf_label() 

#' We use spdep::knearneigh to ensure we obtain 4 of the nearest neighbors to a polygon
#' 
#' However, this system has limitations, if an island is ringed by land, it would work adequately.
#' But when the island is off a coast, it may return neighbors which are 'behind' another neighbor. 
#' We draw lines from the focal portion of the grid we will be combining to each neighbor, 
#' if the lines cross any other neighbor they are discarded. 
#' 
#' This is preferable to simply using the nearest feature (e.g. sf:;st_nearest_feature) when the chain of islands is elongated and may more appropriately be split across multiple grids downstream from here. 
#' @param from a point on surface of the grid which will be merged. POS will ensure it lands on a feature.
#' @param destinations the set of kearneigh which are accepting polygon merges 
first_near_neigh <- function(from, destinations){
  
  sf::st_union(from, destinations[x,]) |> 
    sf::st_cast('LINESTRING') |>
    sf::st_intersects( merge_into_sf[neighbors[[]], ]) |>
    lengths()
  
}








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
    prop_target[[i]][!kp] <-prop_areas[[i]][!kp]
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

gr <- dplyr::mutate(gr, ID = 1:dplyr::n(),  .before = x) |>
  dplyr::rename(geometry = x)
neighb_grid <- vector(mode = 'list', length = length(prop_areas))
for (i in seq_along(neighbors)){
  neighb_grid[[i]] <- gr[neighbors[[i]], ]
}

# place points throughout the grids which need to be merged to determine
# how they will be reallocated into larger grids. 
to_merge_sf <- dplyr::rename(to_merge_sf, geometry = x)
to_merge_sf <- split(to_merge_sf, f = 1:nrow(to_merge_sf))

out <- vector(mode = 'list', length = length(prop_areas))
for (i in seq_along(out)){
  out[[i]] <- assignGrid_pts(
    neighb_grid =  neighb_grid[[i]], 
    focal_grid = to_merge_sf[[i]], 
    props = props[[i]], 
    nf_pct = nf_pct[[i]]
  )
}

# finally create polygons from the point samples
final <- vector(mode = 'list', length = length(out))
for (i in seq_along(final)){
  final[[i]] <- snapGrids(
    x = out[[i]],
    neighb_grid =  neighb_grid[[i]], 
    focal_grid = to_merge_sf[[i]]
  )
} 
final_grids <- dplyr::bind_rows(final)

groups <- split(final_grids, f = final_grids$Assigned)
final_grids <- lapply(groups, snapR) |> bind_rows()
final_grids <- sf::st_make_valid(final_grids)

# reconstitute all original input grids, i.e. those without neighbors, 
# with all reassigned grids. 
gr2 <- sf::st_difference(
  gr,  sf::st_make_valid(sf::st_union(sf::st_combine(final_grids)))
) |>
  sf::st_make_valid()
gr2 <- gr2[as.numeric(sf::st_area(gr2))/1000 > 10,]
final_grids <- gr2 |>
  dplyr::rename(Assigned = ID) |>
  dplyr::bind_rows(final_grids)  |>
  sf::st_as_sf()

# now number the grids in a uniform fashion
cents <- sf::st_point_on_surface(final_grids)
cents <- cents |>
  mutate(
    X = sf::st_coordinates(cents)[,1],
    Y = sf::st_coordinates(cents)[,2]
  ) |>
  dplyr::arrange(-Y, X) |>
  dplyr::mutate(NEWID = 1:dplyr::n()) |>
  sf::st_drop_geometry() |>
  dplyr::select(Assigned, NEWID)

final_grids <- dplyr::left_join(final_grids, cents, by = 'Assigned') |>
  dplyr::select(-Assigned, Assigned = NEWID)

