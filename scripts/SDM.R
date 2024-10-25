# install.packages('dismo', 'spdep', 'spThin', 'glmnet', 'caret', 'CAST')
setwd('~/Documents/assoRted/StrategizingGermplasmCollections')
source('./scripts/createPCNM_fitModel.R')
source('./scripts/WriteSDMresults.R')
source('./scripts/postProcessSDM.R')
################################################################################

x <- read.csv(file.path(system.file(package="dismo"), 'ex', 'bradypus.csv'))
x <- x[,c('lon', 'lat')]
x <- dplyr::distinct(x, .keep_all = )

files <- list.files(
  path = file.path(system.file(package="dismo"), 'ex'), 
  pattern = 'grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the independent variables

# Step 0 define spatial extent of study area. 

pts_plan <-   sf::st_transform(
  sf::st_as_sf(x, coords = c('lon', 'lat'), crs = 4326), 
  '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')

bb <- sf::st_bbox(pts_plan)
buff_dist <- as.numeric( # here we get the mean distance of the XY distances of the bb
  ((bb[3] - bb[1]) + (bb[4] - bb[2])) / 2 
) / 2 # the mean distance * 0.25 is how much we will enlarge the area of analysis. 

bb1 <- sf::st_union(pts_plan) |>
  sf::st_buffer(buff_dist) |>
  terra::vect() |>
  terra::project(terra::crs(predictors)) |>
  terra::ext()

p1 <- terra::mask(predictors, bb1)

rm(pts_plan, bb, buff_dist, bb1)
# Step 1 Select Background points - let's use SDM package envidist for this

pa <- sdm::background(x = p1, n = nrow(x), sp = x, method = 'eDist') |>
  dplyr::select(lon = x,  lat = y)

pa$occurrence <- 0 ; x$occurrence <- 1
x <- dplyr::bind_rows(x, pa) |> # combine the presence and pseudoabsence points
  sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)  |>
  dplyr::mutate(occurrence = factor(occurrence))

brady.df <- data.frame(Species = 'Species', data.frame(sf::st_coordinates(x)))

dists <- sf::st_distance(x[ sf::st_nearest_feature(x), ], x, by_element = TRUE)
thinD <- as.numeric(quantile(dists, c(0.05)) / 1000) # ARGUMENT TO FN @PARAM 

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

index <- unlist(caret::createDataPartition(x$occurrence, p=0.8)) # @ ARGUMENT TO FN @PARAM
train <- x[index,]
test <- x[-index,]

rm(index)

# Develop CV folds for modelling
indices_knndm <- CAST::knndm(train, predictors, k=5)

# Recursive feature elimination using CAST developed folds

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

rm(train1, train_dat)

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

obs <- createPCNM_fitModel(
    x = train, 
    planar_proj = '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')

mod <- obs$mod; cv_model <- obs$cv_model; pcnm <- obs$pcnm

predictors <- c(predictors, pcnm)

rm(train)
# get the variables to extract from the rasters for creating a matrix for 
# predictions, glmnet predict is kind of wonky and needs exact matrix dimensions. 

vars <- rownames(coef(mod)); vars <- vars[2:length(vars)]

# now we need just the COORDINATES FOR TEST and will extract the data from
# this set of predictors to them... 
predict_mat <- predictors[[vars]]
predict_mat <- as.matrix(
  terra::extract(predict_mat, test, ID = FALSE) 
)

cm <- caret::confusionMatrix(
  data = as.factor(predict(mod, newx = predict_mat, type = 'class')), 
  reference = test$occurrence,
  positive="1")

rm(predict_mat)
## Predict our model onto a gridded surface (raster) ## This will allow for downstream
# use with the rest of the safeHavens workflow. 
preds <- predictors[[vars]]
predfun <- function(model, data, ...){
  predict(model, newx=as.matrix(data), type = 'response')
}

rast_cont <- terra::predict(preds, model = mod, fun=predfun, na.rm=TRUE)

rm(lmProfile, predfun, preds)

ob <- postProcessSDM(rast_cont, thresh_metric = 'sensitivity', quant_amt = 0.25)
f_rasts <- ob$f_rasts
thresh <- ob$thresh

terra::plot(f_rasts)
writeSDMresults(
  file.path( 'results', 'SDM'), 'Bradypus_test')



identifyClusters <- function(f_rasts, predictors){
  
  # check if any coefficients are shrunk out of the model, they will be a 0.0000
  # we will remove these from the analysis. 
  
  if(any(unlist(coef(mod))==0)){
    
    retained_terms <- which(as.numeric(coef(mod))!=0)
    retained_terms <- retained_terms[2:length(retained_terms)]
    vars <- rownames(coef(mod)) # all variables to subset from raster stack
    vars <- vars[retained_terms] # remove any shrunk variables
    abs_coef <- abs(c(as.numeric(coef(mod))[retained_terms], 0.1, 0.1)) 
    
  } else {
    
    vars <- rownames(coef(mod)); vars <- vars[2:length(vars)]
    abs_coef <- abs(c(as.numeric(coef(mod)), 0.1, 0.1)) # 
    abs_coef <- abs_coef[2:length(abs_coef)] # remove the intercept term
    
  }
  
  preds <- predictors[[vars]]
  preds$x <- terra::init(preds, fun = 'x') 
  preds$y <- terra::init(preds, fun = 'y') 
  
  pts <- terra::spatSample(
    f_rasts[['Supplemented']], 
    as.points = TRUE,
    method = 'random', size = 500, na.rm = TRUE)
  
  pts <- terra::extract(
    preds, pts,  bind = TRUE
  ) |> 
    as.data.frame() |>
    dplyr::select(-Supplemented)
  
  # add XY and set arbitrary ranks
  
  stanDev <- terra::global(preds, 'sd', na.rm = TRUE)$sd
  abs_coef <- abs_coef * stanDev
  weighted_mat <- sweep(pts, 2, abs_coef, FUN="*")
  
  return(list(
    weighted_mat = weighted_mat, 
    preds = preds, 
    abs_coef = abs_coef, 
    pts = pts))
  
}

predictors <- terra::disagg(predictors, fact = 2)

ic_res <- identifyClusters(f_rasts, predictors = predictors)

weighted_mat <- ic_res$weighted_mat
preds <- ic_res$preds
abs_coef <- ic_res$abs_coef
pts <- ic_res$pts

w_dist <- dist(weighted_mat, method = 'euclidean')
clusters <- hclust(w_dist, method = 'ward.D2')
clusterCut <- cutree(clusters, 20)


# we can also just have nbclust suggest the optimal number of clusters. 
NoClusters <- NbClust::NbClust(
  data = weighted_mat, diss = w_dist, 
  distance = NULL, min.nc = 5, max.nc = 20, 
  method = 'kmean', index = 'silhouette'
)


# Prepare data for training the KNN classifier #
weighted_mat$ID <- factor(clusterCut)

weighted_mat <- weighted_mat[complete.cases(weighted_mat),]
index <- unlist(caret::createDataPartition(weighted_mat$ID, p=0.85)) # @ ARGUMENT TO FN @PARAM
train <- weighted_mat[index,]
test <- weighted_mat[-index,]

# next we will use a split which ensures that a few members of each class 
# are represented in each fold

trainControl <- caret::trainControl(
  method="repeatedcv", number=10, repeats=5)

fit.knn <- caret::train(ID ~ ., data=train, method="knn",
                        trControl = trainControl, metric = 'Accuracy')
knn.k1 <- fit.knn$bestTune # keep this Initial k for testing with knn() function in next section

predicted <- predict(fit.knn, newdata = test)
confusionMatrix(predicted, test$ID)

preds <- preds*abs_coef # need to put onto the rescaled... scale. 
out <- terra::predict(preds, model = fit.knn, na.rm = TRUE)

# outProbs <- terra::predict(preds, model = fit.knn, na.rm = TRUE)
out <- terra::mask(out, f_rasts$Supplemented)
terra::plot(out)

terra::writeRaster(out, './results/SDM/clusterTest.tif', overwrite = TRUE)



################################################################################
############ IF DATA SET IS HEAVILY UNBALANCED THEN DO THIS  ###################

# we surely need to have SOMEWHAT of a balanced data set for the classifier

# let's sample our original raster again, but this time, let's search in the
# geographic spaces where we obtained these groups members from - then we'll
# run the clustering process again and hope we get vaguely similar groups with 
# similar sample sizes. 

# any point with fewer than the median number of observations will be feed back in
cc <- table(clusterCut)
more_samples <- as.numeric(which(cc < median(cc))) # these need more sample


# determine how large each cell is in m, we can use this as a basis for the buffering process. 
r_projected <- f_rasts[['Supplemented']][[1]] |>
  terra::project(
    '+proj=laea +lon_0=-421.171875 +lat_0=-16.8672134 +datum=WGS84 +units=m +no_defs')

d <- terra::xres(r_projected)

# need to use the untransformed points to get the proper sizes. 
need_more_samples <- pts[ weighted_mat$ID %in% more_samples, c('x', 'y')] |>
  sf::st_as_sf(coords = c(x='x', y='y'), crs = terra::crs(f_rasts[['Supplemented']]))  |>
  sf::st_transform(  crs = terra::crs(r_projected)) |>
  sf::st_buffer(d*3.05) |>  # @PARAM IN FUNCTION - HOW MUCH INCREASE AROUND FOCAL CELL!?
  dplyr::summarize(geometry = sf::st_union(geometry)) |>
  sf::st_make_valid() 

concentrated_pts <- sf::st_sample(need_more_samples, size = 100, type = 'regular') |>
  sf::st_as_sf() |>
  sf::st_transform( terra::crs(f_rasts[['Supplemented']]) )

ggplot() +
  geom_sf(data = need_more_samples) +
  geom_sf(data = concentrated_pts)

concentrated_pts <- terra::extract(
  preds, concentrated_pts,  bind = TRUE, 
) |>
  as.data.frame()

concentrated_pts <- concentrated_pts[complete.cases(concentrated_pts),]
concentrated_pts <- unique(concentrated_pts)

weighted_mat2 <- sweep(concentrated_pts, 2, abs_coef, FUN="*")
