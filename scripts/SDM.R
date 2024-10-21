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

rm(thinned, brady.df)

# Step 1.3 - Extract data to points for modelling
bradypus <- terra::extract(predictors, bradypus, bind = TRUE) |>
  sf::st_as_sf() 

# Step 1.2 - create a data split for testing the residuals of the glmnet model
# It's not ideal to do a simple split of these data, because SAC will mean that
# our results could be overly optimistic. SO we don't report these results, 
# we only use them to calculate the residuals from glmnet to then
# calculate MORANS I to determine the effect of spatial
# autocorrelation on the model. 

index <- unlist(caret::createDataPartition(bradypus$occurrence, p=0.8))
train <- bradypus[index,]
test <- sf::st_drop_geometry(bradypus[-index,])

# Fit a simple model to the data and determine whether Spatial autocorrelation is
# present in the residuals. If so, we will apply spThin. 

model <- glm(factor(occurrence) ~ . , data = sf::st_drop_geometry(train), family = 'binomial')
nb <- spdep::knearneigh(train, 4) # now create a neighbor object between the points
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), train, type="idw")
spdep::moran.test(model$residuals, lw)

# Step 2 Develop CV folds for steps 3 and 4

library(NbClust)

chc <- hclust(dist(sf::st_coordinates(train)), method="complete")
train$Cluster <- cutree(chc, k = 5)

indices_LLO <- CAST::CreateSpacetimeFolds(train, spacevar = 'Cluster', k=5)
indices_knndm <- CAST::knndm(train, predictors, k=5)

# Step 3. Recursive feature elimination using CAST developed folds
ctrl <- caret::rfeControl(
  method = "LG0CV",
  repeats = 5,
  number = 10,
  functions = caret::lrFuncs,
  index = indices_knndm$indx_train,
  verbose = FALSE)

lmProfile <- caret::rfe(
  method = 'glmnet',
  sizes = c(2:9), 
  sf::st_drop_geometry(train)[,2:10],
  sf::st_drop_geometry(train)$occurrence,
  rfeControl = ctrl)

predictors(lmProfile) # this is how we subset the relevant variables. 

# Step 4. Model fitting using CAST developed folds
train1 <- dplyr::mutate(
  train, 
  occurrence = dplyr::if_else(occurrence==1, 'YES', 'NO'))

test_class_cv_model <- train(
  x = sf::st_drop_geometry(train1[,predictors(lmProfile)]), 
  st_drop_geometry(train)$occurrence, 
  method = "glmnet", 
  family = 'binomial',
  index = indices_knndm$indx_train)

test_class_cv_model

sub <- sf::st_drop_geometry(train[,predictors(lmProfile)])

# now fit the model just using glmnet::glment in order that we can get the 
# type of response for type='prob' rather than log odds or labelled classes
# which we need to work with terra::predict. 
mod <- glmnet::glmnet(
  x = sub, 
  st_drop_geometry(train)$occurrence, 
  family = 'binomial', 
  keep = TRUE,
  lambda = 0.04821905, alpha = 0.1
)

coef(mod)
# to get the residuals we can do the following

varImp(mod, mod$lambda)
jvars <- as.matrix(test[, predictors(lmProfile)])
predictions_LASSO <- predict(mod, newx = jvars, type = 'class')

confusionMatrix( as.factor(predictions_LASSO), as.factor(test$occurrence))
library(glmnet)
library(ggfortify)
autoplot(mod, colour = 'blue')


nb <- spdep::knearneigh(bradypus, 4) # now create a neighbor object between the points
# the neighbors are scaled by the inverse distance weight between them. 
lw <- spdep::nb2listwdist(spdep::knn2nb(nb), bradypus, type="idw")
spdep::moran.test(mod$residuals, lw) # determine spatial autocorrelation within the model

#preds <- predictors[[c("biome", "bio1", 'bio16', 'bio17')]]
pred <- terra::predict(predictors, mod, type = 'response')
predictors







x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
g2=sample(1:2,100,replace=TRUE)
g4=sample(1:4,100,replace=TRUE)
fit1=glmnet::glmnet(x,y)
predict(fit1,newx=x[1:5,],s=c(0.01,0.005))
predict(fit1,type="coef")
fit2=glmnet::glmnet(x,g2,family="binomial")
predict(fit2,type="response",newx=x[2:5,])
predict(fit2,type="nonzero")
fit3=glmnet(x,g4,family="multinomial")
predict(fit3,newx=x[1:3,],type="response",s=0.01)






library(ggfortify)









library(glmmTMB)
Owls <- transform(Owls,
                  Nest=reorder(Nest,NegPerChick),
                  NCalls=SiblingNegotiation,
                  FT=FoodTreatment)

fit_zipoisson <- glmmTMB(NCalls~(FT+ArrivalTime)*SexParent+
                           offset(log(BroodSize))+(1|Nest),
                         data=Owls,
                         ziformula=~1,
                         family=poisson)



train$pos <- numFactor(train$x, train$y)
# then create a dummy group factor to be used as a random term
train$ID <- factor(rep(1, nrow(train)))

# fit the model
m_tmb <- glmmTMB(occurrence ~ bio17 + bio12 + bio16 + bio1 + exp(pos + 0 | ID), train,
                 family = 'binomial') # take some time to fit

summary(m_tmb)
sims <- DHARMa::simulateResiduals(m_tmb)
plot(sims)


