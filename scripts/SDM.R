library(dismo)
library(terra)
library(spdep)

################################################################################
bradypus <- read.csv(paste0(system.file(package="dismo"), "/ex/bradypus.csv"))
bradypus <- bradypus[,c('lon', 'lat')]

files <- list.files(path=paste(system.file(package="dismo"), '/ex',
                               sep=''),  pattern='grd',  full.names=TRUE )
predictors <- terra::rast(files) # import the indepedent variables
predictors$x <- init(predictors, fun = 'x') 
predictors$y <- init(predictors, fun = 'y') 
# Step 1 Select Background points - let's use SDM package envidist for this

pa <- sdm::background(x = predictors, n = nrow(bradypus), sp = bradypus, method = 'eDist') |>
  dplyr::select(lon = x,  lat = y)

pa$occurrence <- 0 ; bradypus$occurrence <- 1
bradypus <- dplyr::bind_rows(bradypus, pa) |> # combine the presence and pseudoabsence points
  sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)  |>
  dplyr::mutate(occurrence = factor(occurrence))

brady.df <- data.frame(Species = 'Species', data.frame(sf::st_coordinates(bradypus)))

dists <- sf::st_distance(bradypus[ sf::st_nearest_feature(bradypus), ], bradypus, by_element = TRUE)
thinD <- as.numeric(quantile(dists, c(0.1)) / 1000)

thinned <- spThin::thin(loc.data = brady.df, thin.par = thinD,
             spec.col = 'Species',
             lat.col = 'Y', long.col = 'X', reps = 100, 
             locs.thinned.list.return = TRUE, 
             write.files = FALSE, 
             write.log.file = FALSE)

thinned <- data.frame(thinned[ which.max(unlist(lapply(thinned, nrow)))]) |>
  sf::st_as_sf(coords = c('Longitude', 'Latitude'), crs = 4326)

bradypus <- bradypus[lengths(sf::st_intersects(bradypus, thinned))>0,]

rm(thinned, brady.df, pa, thinD, dists, files)

# Step 1.3 - Extract data to points for modelling
bradypus <- terra::extract(predictors, bradypus, bind = TRUE) |>
  sf::st_as_sf() 

# Step 1.2 - create a data split for testing the residuals of the glmnet model
# It's not ideal to do a simple split of these data, because SAC will mean that
# our results could be overly optimistic. SO we don't report these results, 
# we only use them to calculate the residuals from glmnet to then
# calculate MORANS I to determine the effect of spatial
# autocorrelation on the model. 

index <- unlist(caret::createDataPartition(bradypus$occurrence, p=0.85))
train <- bradypus[index,]
test <- sf::st_drop_geometry(bradypus[-index,])

rm(index)
# Fit a simple model to the data and determine whether Spatial autocorrelation is
# present in the residuals. If so, we will apply spThin. 

model <- glm(factor(occurrence) ~ . , data = sf::st_drop_geometry(train), family = 'binomial')
nb <- spdep::knearneigh(train, 4) # now create a neighbor object between the points
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), train, type="idw")
morI <- spdep::moran.test(model$residuals, lw)


predict(model)[1:5] #- predict(model, type = 'response') #== model$residuals
vals <- (as.numeric(train$occurrence) - model$fitted.values) / (model$fitted.values * (1 - model$fitted.values))
options(scipen = 999)

observed <- as.numeric(train$occurrence) - 1
predicted <- predict(model, newdata = train, type = 'response') 
residuals <- observed - predicted

rm(nb, lw, model)
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
  sizes = c(2:ncol(train_dat)), 
  train_dat,
  sf::st_drop_geometry(train)$occurrence,
  rfeControl = ctrl)


rm(ctrl)
# Step 4. Model fitting using CAST developed folds
train1 <- dplyr::mutate(
  train, # requires a YES OR NO or T/F just not anything numeric alike. 
  occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))

cv_model <- train(
  x = sf::st_drop_geometry(train1[,predictors(lmProfile)]), 
  st_drop_geometry(train)$occurrence, 
  method = "glmnet", 
  family = 'binomial',
  index = indices_knndm$indx_train)

sub <- train_dat[,predictors(lmProfile)]

rm(indices_knndm, train1)
# now fit the model just using glmnet::glment in order that we can get the 
# type of response for type='prob' rather than log odds or labelled classes
# which we need to work with terra::predict. 
mod <- glmnet::glmnet(
  x = sub, 
  st_drop_geometry(train)$occurrence, 
  family = 'binomial', 
  keep = TRUE,
  lambda = cv_model$bestTune$lambda, alpha = cv_model$bestTune$alpha
)


rm(cv_model, sub)
# get model information below
coef(mod)
varImp(mod, mod$lambda)
predict_mat <- as.matrix(test[, predictors(lmProfile)])

confusionMatrix(
  as.factor(predict(mod, newx = predict_mat, type = 'class')), 
  as.factor(test$occurrence))

ob <- predict(mod, newx = predict_mat)

nb <- spdep::knearneigh(train, 4) # now create a neighbor object between the points
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), train, type="idw")
# rm(test, train, predict_mat)


preds <- predictors[[predictors(lmProfile)]]
predfun <- function(model, data, ...){
  predict(model, newx=as.matrix(data), type = 'response')
}

x <- terra::predict(preds, model = mod, fun=predfun, na.rm=TRUE)
plot(x)

rm(lmProfile, predfun)






# determine a threshold for creating a binomial map of the species distribution
# we want to predict MORE habitat than exists, so we want to maximize sensitivity
# in our classification. 

test.sf <- sf::st_as_sf(test, coords = c('x', 'y'), crs = 4326) |>
  dplyr::select(occurrence)

test.sf <- terra::extract(x, test.sf, bind = TRUE) |>
  sf::st_as_sf() |>
  sf::st_drop_geometry() 


eval_ob <- dismo::evaluate(
  p = test.sf[test.sf$occurrence==1,'s0'],
  a = test.sf[test.sf$occurrence==0,'s0']
)
thresh <- threshold(eval_ob)
cut <- thresh[['sensitivity']]

m <- matrix(
  c(
    0, cut, 0,
    cut, 1, 1), 
  ncol = 3, byrow= TRUE)

x_binary <- terra::classify(x, m) # create a YES/NO raster

rm(eval_ob, thresh, cut, m)
# use sf::st_buffer() to only keep habitat within XXX distance from known populations
plot(x_binary)
points(bradypus[bradypus$occurrence==1,])

sf::st_buffer()

