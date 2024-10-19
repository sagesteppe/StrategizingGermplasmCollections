setwd('~/Documents/assoRted/StrategizingGermplasmCollections/scripts')

ecoregions <- sf::st_read('../data/spatial/us_eco_l4/us_eco_l4_no_st.shp', quiet = TRUE) |>
  sf::st_transform(4326) |>
  sf::st_make_valid()

polygon <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME == 'Colorado') |>
  sf::st_transform(4326)

EcoregionBasedSample <- function(x, ecoregions, n, increase_method, reduction_method){
  
  if(missing(n)){n<-20} 
  if(missing(increase_method)){increase_method<-'Area'}
  if(missing(reduction_method)){reduction_method<-'Largest'}
  sf::st_agr(ecoregions) = 'constant'
  sf::st_agr(x) = 'constant'
  
  # reduce their input ecoregions data set to the areas where x is. 
  ecoregions <- sf::st_intersection(ecoregions, x) |>
    sf::st_make_valid() |>
    dplyr::mutate(ID = 1:dplyr::n(), .before = 1)
  
  area <- dplyr::mutate(
    ecoregions, Area = units::set_units(sf::st_area(ecoregions), ha), .before = geometry)
  
  area_summaries <- area |>
    sf::st_drop_geometry() |>
    dplyr::group_by(L4_KEY) |>
    dplyr::summarise(
      Eco_lvl = 4,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = L4_KEY) |>
    dplyr::ungroup()
  
  eco_lvls_ct <- data.frame(
    Eco_lvl = c('L1', 'L2', 'L3', 'L4'), 
    ct = c(length(unique(area$L1_KEY)), length(unique(area$L2_KEY)),
           length(unique(area$US_L3CODE)), length(unique(area$US_L4CODE)))
  )
  
  cols <- c('ID', 'US_L4CODE', 'US_L4NAME', 'n', 'geometry')
  
  # in this method we only assign counts to each area. 
  if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] == n){
    
    if(sum(area_summaries[area_summaries, 'Polygon_ct']) == n){
      # if n polygons == n  woohoo! each ecoregion is allocated a sample size of one point 
      
      out <- dplyr::mutate(polygons, n = 1) |>
        dplyr::select(cols)
      
    } else {
      # if n polygons > n, select a polygon in each ecoregion to represent it.  Either by AREA or CENTRALITY 
      if(increase_method == 'Area'){
        
        out <- dplyr::group_by(polygons, US_L4CODE) |>
          dplyr::arrange(Area, .by_group = TRUE) |>
          slice_max(n = 1) |>
          dplyr::mutate(n = 1) |>
          dplyr::select(cols)
        
      } else {
        
        unions <- dplyr::group_by(polygons, US_L4CODE) |>
          dplyr::summarise(geometry = sf::st_union(geometry))
        pts <- sf::st_point_on_surface(unions) 
        out <- polygons[lengths(sf::st_intersects(pts, unions))>0, ] |>
          dplyr::mutate(n = 1) |>
          dplyr::select(dplyr::all_of(cols))
        
      }
    }
    
  } else if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] < n) {
    # each L4 ecoregion is represented by a polygon, fewer than n L4's also means 
    # fewer than n polygons. The only option in this scenario is to add multiple
    # points per polygon based on AREA. 
    
    pct_area <- as.numeric(area$Area / sum(area$Area)) * 100
    sample <- numeric(length = length(pct_area))
    sample[pct_area<5] <- 1 # these small areas by default get a point. 
    n_remain <- n - sum(sample)
    
    if(n_remain >= n/2){
      pct_record <- sum(pct_area[pct_area>5]) / n_remain
      sample[pct_area>5] <- round(pct_area[pct_area>5] / pct_record)
      out <- dplyr::mutate(area, n = sample) |>
        dplyr::select(dplyr::all_of(cols))
    } else {
      # some areas bay be composed entirely of very small coverage areas. 
      # the largest 20 will get the records in that case. 
      area$n <- NA
      area$n[(order(pct_area, decreasing = TRUE)[1:20])] <- 1
      out <- dplyr::select(area, dplyr::all_of(cols))
    }
    
  } else { # MANY L4's across the species range, 
    
    # offer three options for how to select the target L4's
    # 1 & 2, by area, either descending - the largest ones, or ascending - the smallest ones
    # by the number of unique polygons per L4, which can be roughly considered as
    # representing the amount of discontinuity of each L4. 
    
    if(reduction_method=='Largest'){
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Total_area, decreasing = TRUE)[1:n],]$Name, ]
      out <- dplyr::slice_max(out, n = 1, order_by = dplyr::desc(Area), by = US_L4NAME)
    } else if(reduction_method=='Smallest'){
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Total_area, decreasing = TRUE)[1:n],]$Name, ]
      out <- dplyr::slice_max(out, n = 1, order_by = Area, by = US_L4NAME)
    } else {
      out <- area[area$L4_KEY %in% area_summaries[order(area_summaries$Polygon_ct, decreasing = TRUE)[1:n],]$Name, ]
      out <- dplyr::slice_max(out, n = 1, order_by = Area, by = US_L4NAME)
    }
    
    out <- dplyr::mutate(out, n = 1) |>
      dplyr::select(dplyr::all_of(cols))
  }
  
  # above we only added the desired sample sizes to each targeted ecoregion, 
  # now we will add back the ecoregions which are not targeted for sampling to reach
  # the n = 20 goal. 
  
  ids <- unique(dplyr::pull(out, ID))
  out <- ecoregions |>
    dplyr::filter(! US_L4NAME %in% ids) |>
    dplyr::bind_rows(out) |>
    dplyr::arrange(ID) |>
    dplyr::select(dplyr::all_of(cols)) |>
    dplyr::select(-ID) |>
    dplyr::mutate(
      n = dplyr::if_else(is.na(n), 0, n)
    )
  
}


ob <- EcoregionBasedSample(polygon, ecoregions, n = 20)
sum(ob$n)

ggplot() + 
  geom_sf(data = ob, aes(fill = n))
