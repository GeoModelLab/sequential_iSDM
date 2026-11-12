###################################################################################################
############ Algorithms selection (b): Perform cross validation and select best models ############
###################################################################################################

library(terra)
library(blockCV)
library(fuzzySim)
library(modEvA)
library(predicts)
library(maxnet)
library(mgcv)
library(dbarts)
library(randomForest)

### load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
predictorsCropped <- rast(file.path('output', fp_base, 'modellingData/crossValidation/CV_predictorsCropped.tif'))
modData <- read.csv(file.path('output', fp_base, 'modellingData/crossValidation/CV_modData.csv')) # train dataset
modDataSPAT <- vect(modData, geom = c('x', 'y'), crs = 'EPSG:4326') # spatialize it

region <- vect(file.path('output', fp_base, 'modellingData/ModelRegion/CV_modRegion.shp')) # modelling region polygon
world <- geodata::world(path = 'utils/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

vars <- read.csv(file.path('output', fp_base, 'modellingData/selectedVariables.csv'))[, 2]
######### 


######## Blocks Creation --------
## create blocks and define the number of folds for cross-validation
dataSPATsf <- sf::st_as_sf(modData, coords = c('x', 'y'), crs = 'EPSG:4326') # transform to sf object
dataSPATsf <- sf::st_transform(dataSPATsf, crs = "+proj=moll") # use Mollweide projection for euclidean distances
## use spatial autocorrelation of points to select blocks size quantitatively
spatAutocorrPoints <- cv_spatial_autocor(x = dataSPATsf, column = 'presence',
                                         plot = T, progress = T)
blockSize <- spatAutocorrPoints$range
rm(dataSPATsf) # remove sf object because not needed
## blocks area in km
message('The block size in km is: ', round(blockSize/1000, 2), ' km')
## clean unused memory
gc()

## define the k folds for cv and divide the region in blocks of k classes
k = 10
## create the blocks and divide them in five folds
blocks <- cv_spatial(modDataSPAT,
                     column = 'presence', r = predictorsCropped, size = blockSize, 
                     k = k, selection = 'random', hexagon = T, seed = 123)
## assign the belonging fold to each observation in the data (also the spatial one for plotting)
modData$fold <- blocks$folds_ids
modDataSPAT$fold <- blocks$folds_ids

## check the number of test data for each fold
(n <- table(modData[, c('fold', 'presence')]))
(nt <- rowSums(n))  # n total rows per fold
(np <- n[ , 2])  # n presence rows per fold

## plot the folds and presence points
plot(modDataSPAT, main = 'block cross validation', font.main = 1,
     y = 'fold', cex = 0.4, col = hcl.colors(k, 'viridis'), legend = F,
     mar = c(5.1, 2.1, 4.1, 4.1))
plot(world, add = T)
plot(vect(blocks$blocks), border = 'blue', add = TRUE, lw = .3)
legend(legend = sort(unique(modDataSPAT$fold)), col = hcl.colors(k, 'viridis'),
       'right', xpd = T, bty = 'n', inset = c(-0.1, 0),
       title = 'folds', pch = 16, pt.cex = 1, y.intersp = 0.5, x.intersp = 0.2)
mtext(paste('block size:', round(blockSize/1000, 2), 'km'),
     side = 1, line = -1)
## save the blocked data to csv
write.csv(data.frame(modData), file.path('output', fp_base, 'modellingData/crossValidation/CV_foldedData.csv'), row.names = F)

## clean unused memory
gc()
######### 


######## Block Cross Validation training --------
foldedData <- read.csv(file.path('output', fp_base, 'modellingData/crossValidation/CV_foldedData.csv'))
folds <- sort(unique(foldedData$fold))
## create a df to store cv prediction results
predsCV <- foldedData[ , c(1:4, grep("fold", names(foldedData)))]

## iterate through the folds and train the models
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')
for (f in folds) {
  startTime <- proc.time() # start time of the fold
  ## train the model on all data except fold f, but do the prediction using all data
  foldTrainingData <- subset(foldedData, foldedData$fold != f)
  
  message(paste('fold', f, ' - Bioclim'))
  bcMod_fold <- envelope(x = subset(foldTrainingData, presence == 1, select = vars))
  bcPred_fold <- predict(bcMod_fold, foldedData)
  predsCV[ , paste0('BC_fold_', f, '_pred')] <- as.vector(bcPred_fold)
  predsCV[ , paste0('BC_fold_', f, '_pred_reclass')] <- quantReclass(bcPred_fold)
  
  message(paste('fold', f, ' - Maxent'))
  mxtFormula <- maxnet.formula(foldTrainingData[, 'presence'], foldTrainingData[, vars])
  mxtMod_fold <- maxnet(p = foldTrainingData[, 'presence'], data = foldTrainingData[, vars],
                        f = mxtFormula)
  mxtPred_fold <- predict(mxtMod_fold, foldedData, type = 'cloglog')
  predsCV[ , paste0('MXT_fold_', f, '_pred')] <- as.vector(mxtPred_fold)
  
  message(paste('fold', f, ' - GLM'))
  glmFormula <- reformulate(termlabels = vars, response = 'presence')
  glmMod_fold <- glm(formula = glmFormula, family = binomial, data = foldTrainingData)
  glmPred_fold <- predict(glmMod_fold, foldedData, type = 'response')
  glmFav_fold <- Fav(pred = glmPred_fold, sample.preval = prevalence(model = glmMod_fold))
  predsCV[ , paste0('GLM_fold_', f, '_pred')] <- as.vector(glmPred_fold)
  predsCV[ , paste0('GLM_fold_', f, '_fav')] <- as.vector(glmFav_fold)
  
  message(paste('fold', f, ' - GAM'))
  gamFormula <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = "cs", k = 20)', collapse = '+')))
  gamMod_fold <- gam(gamFormula, family = binomial, data = foldTrainingData, method = 'REML')
  gam.check(gamMod_fold)
  gamPred_fold <- predict(gamMod_fold, foldedData, type = 'response')
  gamFav_fold <- Fav(pred = gamPred_fold, sample.preval = prevalence(model = gamMod_fold))
  predsCV[ , paste0('GAM_fold_', f, '_pred')] <- as.vector(gamPred_fold)
  predsCV[ , paste0('GAM_fold_', f, '_fav')] <- as.vector(gamFav_fold)
  
  message(paste('fold', f, ' - BART'))
  bartMod_fold <- bart(x.train = foldTrainingData[, vars], y.train = foldTrainingData[, 'presence'],
                       keeptrees = TRUE, verbose = FALSE, seed = 1312)
  bartPred_fold <- predict_bart_df(bartMod_fold, foldedData)
  bartFav_fold <- Fav(pred = bartPred_fold, sample.preval = prevalence(model = bartMod_fold))
  predsCV[ , paste0('BART_fold_', f, '_pred')] <- as.vector(bartPred_fold)
  predsCV[ , paste0('BART_fold_', f, '_fav')] <- as.vector(bartFav_fold)
  
  message(paste('fold', f, ' - RF'))
  rfFormula <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod_fold <- randomForest(rfFormula, foldTrainingData, na.action = na.exclude)
  rfPred_fold <- predict(rfMod_fold, foldedData)
  rfPred_fold <- pmin(pmax(rfPred_fold, 0), 1) # clamp to 0,1
  rfFav_fold <- Fav(pred = rfPred_fold, sample.preval = prevalence(model = rfMod_fold))
  predsCV[ , paste0('RF_fold_', f, '_pred')] <- as.vector(rfPred_fold)
  predsCV[ , paste0('RF_fold_', f, '_fav')] <- as.vector(rfFav_fold)
  
  endTime <- proc.time() # end time of the fold
  message('Elapsed fold ', f, ' time: ', (endTime-startTime)[3], ' sec') # elapsed fold time
  message('')
  gc()
}

## save outputs to disk
write.csv(predsCV, file.path('output', fp_base, 'models/crossValidation/predsCV.csv'), row.names = FALSE)
######### 

## offload memory
rm(list = ls())
gc()
