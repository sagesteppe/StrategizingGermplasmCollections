# install.packages('dismo', 'spdep', 'spThin', 'glmnet', 'caret', 'CAST')
################################################################################

x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- dplyr::distinct(x, .keep_all = )

files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables

# Step 1 Select Background points - let's use SDM package envidist for this

pa <- sdm::background(x = predictors, n = nrow(x), sp = x, method = 'eDist') |>
  dplyr::select(lon = x,  lat = y)

pa$occurrence <- 0 ; x$occurrence <- 1
x <- dplyr::bind_rows(x, pa) |> # combine the presence and pseudoabsence points
  sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)  |>
  dplyr::mutate(occurrence = factor(occurrence))

brady.df <- data.frame(Species = 'Species', data.frame(sf::st_coordinates(x)))

dists <- sf::st_distance(x[ sf::st_nearest_feature(x), ], x, by_element = TRUE)
thinD <- as.numeric(quantile(dists, c(0.1)) / 1000) # ARGUMENT TO FN @PARAM 

thinned <- spThin::thin(
  loc.data = brady.df, thin.par = thinD,
  spec.col = 'Species',
  lat.col = 'Y', long.col = 'X', reps = 100, 
  locs.thinned.list.return = TRUE, 
  write.files = FALSE, 
  write.log.file = FALSE)

thinned <- data.frame(thinned[ which.max(unlist(lapply(thinned, nrow)))]) |>
  sf::st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)

x <- x[lengths(sf::st_intersects(x, thinned))>0,]

rm(thinned, brady.df, pa, thinD, dists, files)

# Step 1.3 - Extract data to points for modelling
x <- terra::extract(predictors, x, bind = TRUE) |>
  sf::st_as_sf() 

# Step 1.2 - create a data split for testing the residuals of the glmnet model
# It's not ideal to do a simple split of these data, because SAC will mean that
# our results could be overly optimistic. SO we don't report these results, 
# we only use them to calculate the residuals from glmnet to then
# calculate MORANS I to determine the effect of spatial
# autocorrelation on the model. 

index <- unlist(caret::createDataPartition(x$occurrence, p=0.85)) # @ ARGUMENT TO FN @PARAM
train <- x[index,]
test <- x[-index,]

rm(index)
# Fit a simple model to the data and determine whether Spatial autocorrelation is
# present in the residuals. If so, we will apply spThin. 

model <- glm(factor(occurrence) ~ . , data = sf::st_drop_geometry(train), family = 'binomial')
nb <- spdep::knearneigh(train, 4) # now create a neighbor object between the points
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), train, type="idw")
morI <- spdep::moran.test(model$residuals, lw)

predicted <- predict(model, newdata = train, type = 'response') 
residuals <- (as.numeric(train$occurrence) - 1) - predicted

rm(nb, lw, model, predicted, residuals)
# Step 2 Develop CV folds for steps 3 and 4
indices_knndm <- CAST::knndm(train, predictors, k=5)

# Step 3. Recursive feature elimination using CAST developed folds

train_dat <- sf::st_drop_geometry(train[, -which(names(train) %in% c("occurrence"))])
ctrl <- caret::rfeControl(
  method = "LG0CV",
  repeats = 5,
  number = 10,
  functions = caret::lrFuncs,
  index = indices_knndm$indx_train,
  verbose = FALSE)

lmProfile <- caret::rfe(
  method = 'glmnet',
  sizes = c(3:ncol(train_dat)), 
  x = train_dat,
  y = sf::st_drop_geometry(train)$occurrence,
  rfeControl = ctrl)

# Step 4. Model fitting using CAST developed folds
train1 <- dplyr::mutate(
  train, # requires a YES OR NO or T/F just not anything numeric alike. 
  occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))

cv_model <- train(
  x = sf::st_drop_geometry(train1[,predictors(lmProfile)]), 
  sf::st_drop_geometry(train)$occurrence, 
  method = "glmnet", 
  family = 'binomial', 
  index = indices_knndm$indx_train) 

sub <- train_dat[,predictors(lmProfile)]

# rm(train1, train_dat)

# now fit the model just using glmnet::glment in order that we can get the 
# type of response for type='prob' rather than log odds or labelled classes
# which we need to work with terra::predict. 
mod <- glmnet::glmnet(
  x = sub, 
  sf::st_drop_geometry(train)$occurrence, 
  family = 'binomial', 
  keep = TRUE,
  lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
)

rm(cv_model)
###################

train_planar <- sf::st_transform(
  train, '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs') 

dis <- sf::st_distance(train_planar)
dis <- apply(dis, 2, as.numeric)
xypcnm <- vegan::pcnm(dis)
xypcnm.df <- data.frame(xypcnm$vectors)[,1:20]

pcnmProfile <- caret::rfe( 
  method = 'glmnet', # https://doi.org/10.1111/gean.12054
  x = xypcnm.df, 
  sizes = 1:5,
  sf::st_drop_geometry(train)$occurrence,
  rfeControl = ctrl, 
  index = indices_knndm$indx_train
  )


xypcnm.df <- xypcnm.df[,caret::predictors(pcnmProfile)]
preds <- cbind(sub, xypcnm.df)

if(is.numeric(xypcnm.df)){
  colnames(preds)[length(preds)] <- caret::predictors(pcnmProfile)}

rm(xypcnm, dis, train_planar)
# trying to refit the glmnet

cv_model <- train( # let's extract the model performance for the top alpha/lambda info here. 
  x = preds, 
  sf::st_drop_geometry(train)$occurrence, 
  method = "glmnet", 
  family = 'binomial', 
  index = indices_knndm$indx_train) 

mod <- glmnet::glmnet(
  x = preds, 
  sf::st_drop_geometry(train)$occurrence, 
  family = 'binomial', 
  keep = TRUE,
  lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
)

# to predict onto the confusion matrix, we now need to add the PCNM/MEM values
# for the relevant layers onto our independent test data, this will require us
# to create PCNM raster surfaces

xypcnm.sf <- cbind(xypcnm.df, dplyr::select(train, geometry)) |> 
  sf::st_as_sf()


if(is.data.frame(xypcnm.df)){
  pcnm2raster <- function(x){
    
    fit <- fields::Tps(sf::st_coordinates(xypcnm.sf), x)
    p <- terra::rast(predictors[[1]])
    pcnm <- terra::interpolate(p, fit)
    pcnm <- terra::mask(pcnm, predictors[[1]])
    
    return(pcnm)
  }
  
  pcnm <- lapply(xypcnm.df, pcnm2raster)
  pcnm <- terra::rast(pcnm)
  names(pcnm) <- caret::predictors(pcnmProfile)

} else if(is.numeric(xypcnm.df)){
  
  fit <- fields::Tps(sf::st_coordinates(xypcnm.sf), xypcnm.df)
  p <- terra::rast(predictors[[1]])
  pcnm <- terra::interpolate(p, fit)
  pcnm <- terra::mask(pcnm, predictors[[1]])
  names(pcnm) <- caret::predictors(pcnmProfile)
  
  rm(fit, p)
}  

predictors <- c(predictors, pcnm)
terra::plot(predictors)

rm(xypcnm.sf, pcnm2raster, xypcnm.df, pcnmProfile)

# get model information below
vars <- rownames(coef(mod)); vars <- vars[2:length(vars)]

# now we need just the COORDINATES FOR TEST and will extract the data from
# this set of predictors to them... 
predict_mat <- predictors[[vars]]
predict_mat <- as.matrix(
  terra::extract(predict_mat, test, ID = FALSE) 
)

cm <- confusionMatrix(
  as.factor(predict(mod, newx = predict_mat, type = 'class')), 
  as.factor(test$occurrence))

# determine whether there is strong evidence for spatial autocorrelation in the
# residuals. 
fitted_values <- predict(mod, newx = as.matrix(preds))
residuals <- (as.numeric( sf::st_drop_geometry(train$occurrence))-1) - fitted_values

nb <- spdep::knearneigh(sf::st_as_sf(train, coords = c('x', 'y'), crs = 4326) , 4) 
# now create a neighbor object between the points
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), train, type="idw")
morI <- spdep::moran.test(residuals, lw)

rm(predict_mat, residuals, nb, lw, sub, fitted_values)
## Predict our model onto a gridded surface (raster) ## This will allow for downstream
# use with the rest of the safeHavens workflow. 

preds <- predictors[[vars]]
predfun <- function(model, data, ...){
  predict(model, newx=as.matrix(data), type = 'response')
}

rast_cont <- terra::predict(preds, model = mod, fun=predfun, na.rm=TRUE)

rm(lmProfile, predfun, preds)

# determine a threshold for creating a binomial map of the species distribution
# we want to predict MORE habitat than exists, so we want to maximize sensitivity
# in our classification. 

test.sf <- sf::st_as_sf(test, coords = c('x', 'y'), crs = 4326) |>
  dplyr::select(occurrence)

test.sf <- terra::extract(rast_cont, test.sf, bind = TRUE) |>
  sf::st_as_sf() |>
  sf::st_drop_geometry() 

eval_ob <- dismo::evaluate(
  p = test.sf[test.sf$occurrence==1,'s0'],
  a = test.sf[test.sf$occurrence==0,'s0']
)
thresh <- dismo::threshold(eval_ob)
cut <- thresh[['sensitivity']] # ARGUMENT TO FN @PARAM 

m <- matrix( # use this to reclassiy data to a binary raster
  c( # but more simply, turn the absences into NA for easier masking later on? 
    0, cut, NA,
    cut, 1, 1), 
  ncol = 3, byrow= TRUE)

rast_binary <- terra::classify(rast_cont, m) # create a YES/NO raster

terra::plot(rast_binary)
terra::points(x[x$occurrence==1,])

rm(eval_ob, cut, m, test, test.sf)

# use sf::st_buffer() to only keep habitat within XXX distance from known populations
# we'll use another set of cv-folds based on all occurrence data 
# Essentially, we will see how far the nearest neighbor is from each point in each
# fold

nn_distribution <- function(x, y){
  ob <- unlist(x)
  
  nf <- sf::st_distance(
    y[sf::st_nearest_feature(y[ob, ]), ],
    y[ob, ], by_element = TRUE
  )
}

pres <- x[ x$occurrence==1, ]
pres <- sf::st_transform(
  pres, 
  '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')
indices_knndm <- CAST::knndm(pres, predictors, k=10)

nn_dist <- lapply(indices_knndm[['indx_train']], nn_distribution, y = pres)
dists <- unlist(list(lapply(nn_dist, quantile, 0.25))) # ARGUMENT TO FN @PARAM 

within_dist <- sf::st_buffer(pres, median(dists)) |>
  dplyr::summarize(geometry = sf::st_union(geometry)) |>
  sf::st_simplify() |>
  sf::st_transform(terra::crs(rast_binary)) |>
  terra::vect()

rast_clipped <- terra::mask(rast_binary, within_dist)

rm(train_pres, nn_dist, indices_knndm, nn_distribution, within_dist)
####### IF WE HAVE POINTS WHICH ARE FLOATING IN SPACE - I.E. POINTS W/O  
# SUITABLE HABITAT MARKED, THEN LET'S ADD the same amount of suitable habitat 
# to each of them, that was used as the buffer for clipping suitable habitat to the
# points above. 

pres <- sf::st_transform(pres, terra::crs(rast_binary))
outside_binary <- terra::extract(rast_binary, pres, bind = TRUE) |>
  sf::st_as_sf() |>
  dplyr::filter(is.na(s0)) |>
  sf::st_transform(
    '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs') |>
  sf::st_buffer(min(dists)) |>
  dplyr::summarize(geometry = sf::st_union(geometry)) |>
  dplyr::mutate(occurrence = 1) |>
  terra::vect() |>
  terra::project(terra::crs(rast_binary)) |>
  terra::rasterize( rast_binary, field = 'occurrence')

rast_clipped_supplemented <- max(rast_clipped, outside_binary, na.rm = TRUE)


##########   COMBINE ALL RASTERS TOGETHER FOR A FINAL PRODUCT      #############
f_rasts <- c(rast_cont, rast_binary, rast_clipped, rast_clipped_supplemented)
names(f_rasts) <- c('Predictions', 'Threshold', 'Clipped', 'Supplemented')
terra::plot(f_rasts)

rm(outside_binary, pres, rast_clipped_supplemented, rast_clipped, rast_binary)

ctrl
cv_model

##### WE NEED TO SAVE A WIDE NUMBER OF OBJECTS   ####

# 1) THE MODEL
mod # saveRDS()

# 2) MODEL COEFFICIENTS
coef(mod)

# 3) ALL RASTERS
rast_cont

# 4) evaluation statistics - from normal old school split. 
cm

# 5) THRESHOLD VALUES. 
thresh

# 6) THINNED DATA FOR TRAINING MODEL / CV FOLDS
str(cv_model)


# 8) PCNM layers. 



rm(thresh, mod, rast_cont, rast_binary, rast_clipped)




