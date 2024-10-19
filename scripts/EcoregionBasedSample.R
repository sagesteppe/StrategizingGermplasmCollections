ecoregions <- sf::st_read('../data/spatial/us_eco_l4/us_eco_l4_no_st.shp', quiet = TRUE) |>
  sf::st_transform(4326) |>
  sf::st_make_valid()

nc <- spData::us_states |>
  dplyr::select(NAME) |>
  dplyr::filter(NAME == 'Michigan') |>
  sf::st_transform(4326)

area <- sf::st_intersection(nc, ecoregions)
area <- dplyr::mutate(
  area, Area = units::set_units(sf::st_area(area), ha), .before = geometry)
# Make table of the values required for any downstream steps. 

summarizer <- function(x){
  l1_summary <- area |>
    sf::st_drop_geometry() |>
    dplyr::group_by(L1_KEY) |>
    dplyr::summarise(
      Eco_lvl = 1,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = L1_KEY) |>
    dplyr::ungroup()
  
  l2_summary <- area |>
    sf::st_drop_geometry() |>
    dplyr::group_by(L2_KEY) |>
    dplyr::summarise(
      Eco_lvl = 2,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = L2_KEY) |>
    dplyr::ungroup()
  
  l3_summary <- area |>
    sf::st_drop_geometry() |>
    dplyr::group_by(L3_KEY) |>
    dplyr::summarise(
      Eco_lvl = 3,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = L3_KEY) |>
    dplyr::ungroup()
  
  l4_summary <- area |>
    sf::st_drop_geometry() |>
    dplyr::group_by(L4_KEY) |>
    dplyr::summarise(
      Eco_lvl = 4,
      Polygon_ct = dplyr::n(), 
      Total_area = sum(Area)
    ) |>
    dplyr::rename(Name = L4_KEY) |>
    dplyr::ungroup()
  
  bind_rows(l1_summary, l2_summary, l3_summary, l4_summary)
  
}

area_summaries <- summarizer(area)

eco_lvls_ct <- data.frame(
  Eco_lvl = c('L1', 'L2', 'L3', 'L4'), 
  ct = c(length(unique(area$L1_KEY)), length(unique(area$L2_KEY)),
                length(unique(area$US_L3CODE)), length(unique(area$US_L4CODE)))
)

cols <- c('US_L4CODE', 'US_L4NAME', 'n', 'geometry')

if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] == 20){
  
  if(sum(area_summaries[area_summaries$Eco_lvl==4, 'Polygon_ct']) == 20){
    # if n polygons == 20  woohoo! each ecoregion is allocated a sample size of one point 
    
    out <- dplyr::mutate(
      n = 1
    ) |>
      dplyr::select(cols)
    
  } else {
    # if n polygons > 20, select A polygon in each ecoregion to represent it.  Either by AREA or CENTRALITY 
    if(method=='Area'){
      
    } else {
      
    }
  }
  
} else if(eco_lvls_ct[eco_lvls_ct$Eco_lvl=='L4','ct'] < 20) {
  # each L4 ecoregion is represented by a polygon, fewer than 20 L4's also means 
  # fewer than 20 polygons. The only option in this scenario is to add multiple
  # points per polygon based on AREA. 
  
 
  
  
  } else {
  
  # if n L4 ecoregions > 20; try to maximize the distribution of points across the ecoregions. 
    # we will evaluate this spread using a JOIN COUNT? to reduce spatial autocorrelation 
    # the lower the SA, the most equal the sampling schema. 
  
}


EcoregionBasedSample <- function(x, n, ecoregions){
  
  
  ecoregions <- sf::st_intersection(ecoregions, x)
  
  
  
}
