library(dismo)
library(terra)
library(spdep)
library(spatialreg)

bradypus <- read.csv(paste0(system.file(package="dismo"), "/ex/bradypus.csv"))
bradypus <- bradypus[,c('lon', 'lat')]

files <- list.files(path=paste(system.file(package="dismo"), '/ex',
                               sep=''),  pattern='grd',  full.names=TRUE )
mask <- raster(files[1]) # set an extent to which the background points will be restricted
bg <- setNames(data.frame(dismo::randomPoints(mask, nrow(bradypus))), c('lon', 'lat'))
bg$occurrence <- 0 ; bradypus$occurrence <- 1

bradypus <- dplyr::bind_rows(bradypus, bg) |> # combine the presence and pseudoabsence points
  sf::st_as_sf(coords = c('lon', 'lat'), crs = 4326)  |>
  dplyr::mutate(occurrence = factor(occurrence))

predictors <- terra::rast(files) # import the indepedent variables
bradypus <- terra::extract(predictors, bradypus, bind = TRUE) |>
  sf::st_as_sf() # extract the indepedent variables to the dependent variables

bp.df <- sf::st_drop_geometry(bradypus) # GLM will be confused by an SF object, simplify to data.frame

rm(files, bg, mask) # the environment is too cluttered. 

# create a model where presence is predicted by a few variables. 
glmbase <- glm(occurrence ~ bio1 + bio12 + bio17 + bio5 + biome,  data = bp.df, family=binomial)

nb <- spdep::knearneigh(bradypus, 4) # now create a neighbor object between the points
lw <- spdep::nb2listw(spdep::knn2nb(nb))

spdep::moran.test(glmbase$residuals, lw) # determine spatial autocorrelation within the model

bp.df$preds <- glmbase$fitted.values

MEbinom <- spatialreg::ME(occurrence ~ bio1 + bio12 + bio17 + bio5 + biome,  data = bp.df, family="binomial",
                listw = lw, alpha=0.1, verbose=TRUE, nsim=49) # make and determine 

glmME <- glm(occurrence ~ bio1 + bio12 + bio17 + bio5 + biome + fitted(MEbinom), data = bp.df, family="binomial")

anova(glmbase, glmME, test = 'Chisq')
summary(glmbase)
summary(glmME)

pred <- terra::predict(predictors, glmbase, type = 'response')
plot(pred)


# Step 1 Select Background points - let's use SDM package envidist for this

# Step 2 Develop CV folds for steps 3 and 4

library(NbClust)

coords <- st_coordinates(bradypus)

chc <- hclust(dist(sf::st_coordinates(bradypus)), method="complete")
bradypus$Cluster <- cutree(chc, k = 5)

indices_LLO <- CAST::CreateSpacetimeFolds(bradypus, spacevar = 'Cluster', k=5)
indices_knndm <- knndm(bradypus, predictors, k=5)

plot(geodist(bradypus, predictors, cvfolds = indices_knndm$indx_test))


# Step 3. Recursive feature elimination using CAST developed folds
ctrl <- rfeControl(
  method = "LG0CV",
  repeats = 5,
  number = 10,
  functions = lrFuncs,
  index = indices_knndm$indx_train,
  verbose = FALSE)

lmProfile <- rfe(
  method = 'glmnet',
  sizes = c(1:9),
  st_drop_geometry(bradypus)[,2:10],
  st_drop_geometry(bradypus)$occurrence,
  rfeControl = ctrl)

predictors(lmProfile) # this is how we subset the relevant variables. 

# Step 4. Model fitting using CAST developed folds
bradypus <- dplyr::mutate(
  bradypus, 
  occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))

test_class_cv_model <- train(
  x = sf::st_drop_geometry(bradypus[,predictors(lmProfile)]), 
  st_drop_geometry(bradypus)$occurrence, 
  method = "glmnet", 
  family = 'binomial',
  index = indices_knndm$indx_train)

test_class_cv_model

# now fit the model just using glmnet::glment in order that we can get the 
# type of response for type='prob' rather than log odds or labelled classes
# which we need to work with terra::predict. 

pred <- terra::predict(predictors, test_class_cv_model, type = 'response')

