########################################################################
############ Run the models yearly observations iteratively ############
########################################################################

library(tidyverse)
library(terra)
library(fuzzySim)
library(mgcv)
library(dbarts)
library(embarcadero)
library(randomForest)
library(vip)
library(pdp)
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')

## load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
world <- geodata::world(path = 'utils/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

native <- read.csv(file.path('output', fp_base, 'modellingData/native_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'n')
invasive <- read.csv(file.path('output', fp_base, 'modellingData/invasive_occs.csv')) |> 
  mutate(year = year(eventDate),
         status = 'i') |> 
  ## filter obs before 1996 (first invasive discovery)
  filter(year >= 1996)

predictors <- rast(file.path('output', fp_base,'modellingData/predictors.tif'))
bias <- geodata::footprint(year = 2009, path = 'utils')
bias <- project(bias, predictors)

selectedMods <- read.csv(file.path('output', fp_base, 'models/crossValidation/selectedModels.csv'))
vars <- read.csv(file.path('output', fp_base, 'modellingData/selectedVariables.csv'))[, 2] # variables after collinearity analysis
## modelling region multiplier
wd <- readRDS(file.path('output', fp_base, 'modellingData/ModelRegion/region_mult.rds'))
######### 


######### Loop over the years and run the models ----------
## params
mult_p <- 10                                    # ratio of absences
modData <- list()                               # save all modelling data in a list
allModels <- NULL                               # all models results in raster format
allFavAVG <- NULL                               # all yearly mean predictions
ens <- NULL                                     # initialise ensemble object
regions <- NULL                                 # training regions of all models
allVarimp <- data.frame()                       # variable importance of all models
allRc <- data.frame()                           # response curves

years <- unique(invasive$year)
(years <- sort(years))

# ## plot the observations across all years
# for (y in years) {
#   plot(world, main = paste0('obs in year ', y))
#   points(vect(native, geom = c('long', 'lat'), crs = crs(world)), col = 'blue', cex = 0.2)
#   points(vect(invasive |> filter(year <= y),
#               geom = c('long', 'lat'), crs = crs(world)), col = 'red')
# }

## overall time tracking
(start <- Sys.time())
for (i in 1:length(years)) {
  y <- years[i]
  cat('Iteration', i, ' -  Year', y, '\n')
  ## yearly time tracking
  startTime <- proc.time()
  
  ## modelling data ----
  trainDF <- rbind(invasive[invasive$year <= y, ], native)
  occs <- vect(trainDF, geom = c('long', 'lat'), crs = crs(predictors))

  ## compute the buffer
  b <- buffer(occs, width = wd)
  ## get the region and assign an identifier
  region <- aggregate(b)
  region$year <- y
  ## save the region to output stack
  if (is.null(regions)) {
    regions <- region
  } else {
    regions <- rbind(regions, region)
  }
  ## crop predictors
  predictorsCropped <- crop(predictors, region, mask = T, touches = F)

  ## data gridding
  gridData <- gridRecords(rst = predictorsCropped, pres.coords = occs)
  ## clean from NAs in predictors if any
  gridData <- drop_na(gridData)
  ## crop bias layer for absence generation to study area
  biasCropped <- crop(bias, region, mask = T)
  ## absence generation
  data <- selectAbsences(gridData,
                         sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                         mult.p = mult_p, seed = 1312, plot = F, df = T, bias = biasCropped, verbosity = 0)
  ## save modelling data
  modData[[as.character(y)]] <- data
  #####

  ## model runs ----
  cat('\tGAM ')
  gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = \'cr\')', collapse = '+')))
  gamMod <- bam(gamForm, family = binomial, data = data, model = T,
                discrete = T, nthreads = parallel::detectCores() - 1)
  gamPred <- predict(predictors, gamMod, type = 'response', na.rm = T)
  names(gamPred) <- paste0('gamPred_', y)
  gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
  names(gamFav) <- paste0('gamFav_', y)
  gc()

  cat('- BART ')
  bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
                  ntree = 200, keeptrees = T, verbose = F, seed = 1312)
  bartPred <- predict2.bart(bartMod, raster::stack(predictors), quantiles = c(), splitby = 10, quiet = T)
  bartPred <- rast(bartPred)
  crs(bartPred) <- crs(predictors)
  names(bartPred) <- paste0('bartPred_', y)
  bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
  names(bartFav) <- paste0('bartFav_', y)
  gc()

  cat('- RF\n')
  rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod <- randomForest(rfForm, data, na.action = na.exclude)
  rfPred <- predict(predictors, rfMod, na.rm = T, cores = parallel::detectCores() - 1, cpkgs = 'randomForest')
  names(rfPred) <- paste0('rfPred_', y)
  rfPred <- clamp(rfPred, lower = 0, upper = 1) # clamp to 0,1
  rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
  names(rfFav) <- paste0('rfFav_', y)
  gc()
  
  ## save the models together
  allYear <- c(gamPred, gamFav,
               bartPred, bartFav,
               rfPred, rfFav)
  if (is.null(allModels)) {
    allModels <- allYear
  } else {
    allModels <- c(allModels, allYear)
  }
  #####
  
  
  ## subsample of data & wrapper functions (for vip and pdp) ----
  set.seed(1312)
  subdata <- data |>
    dplyr::slice_sample(n = 5000)
  
  gam_pred_wrapper <- function(object, newdata) {
    p <- predict(object, newdata, type = 'response')
  }
  rf_pred_wrapper <- function(object, newdata) {
    p <- predict(object, newdata)
  }
  bart_pred_wrapper <- function(object, newdata) {
    p <- predict_bart_df(object, newdata)
  }
  #####
  
  ## variable importance ----
  ## gam
  GAMvarimp <- vi_permute(gamMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                          pred_wrapper = gam_pred_wrapper, nsim = 10,
                          parallel = T, .packages = 'mgcv')
  GAMvarimp <- cbind(mod = 'GAM', year = y, GAMvarimp)
  ## rf
  RFvarimp <- vi_permute(rfMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                         pred_wrapper = rf_pred_wrapper, nsim = 10,
                         parallel = T, .packages = 'randomForest')
  RFvarimp <- cbind(mod = 'RF', year = y, RFvarimp)
  ## bart
  BARTvarimp <- vi_permute(bartMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                           pred_wrapper = bart_pred_wrapper, nsim = 10,
                           parallel = T, .packages = 'dbarts')
  BARTvarimp <- cbind(mod = 'BART', year = y, BARTvarimp)
  ## append to results
  allVarimp <- rbind(allVarimp, GAMvarimp, RFvarimp, BARTvarimp)
  cat('\tvariable importance done\n')
  #####
  
  ## response curves ----
  ## gam
  GAMrc <- lapply(vars, function(v) {
    p <- pdp::partial(gamMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, pred.fun = gam_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'GAM',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## rf
  RFrc <- lapply(vars, function(v) {
    p <- pdp::partial(rfMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, prob = T, pred.fun = rf_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'RF',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## bart
  BARTrc <- lapply(vars, function(v) {
    p <- pdp::partial(bartMod, pred.var = v, train = subdata[, vars],
                      grid.resolution = 50, prob = T, pred.fun = bart_pred_wrapper)
    p <- p |> 
      group_by(pick(1)) |> 
      summarise(yhat = mean(yhat),
                .groups = 'drop')
    ## out
    data.frame(mod = 'BART',
               year = y,
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## append to results df
  # allRc <- rbind(allRc, MXTrc, GAMrc, RFrc, BARTrc)
  allRc <- rbind(allRc, GAMrc, RFrc, BARTrc)
  cat('\tresponse curve done\n')
  #####
  
  ## models average ----
  favMean <- app(c(bartFav[[1]], gamFav, rfFav), mean)
  names(favMean) <- paste0('mods_mean_fav_', y)
  favVar <- app(c(bartFav[[1]], gamFav, rfFav), var)
  names(favVar) <- paste0('mods_var_fav_', y)
  favAVG <- c(favMean, favVar)
  if (is.null(allFavAVG)) {
    allFavAVG <- favAVG
  } else {
    allFavAVG <- c(allFavAVG, favAVG)
  }
  #####
  
  ## yearly time tracking
  endTime <- proc.time()
  cat('Elapsed year', y, 'time:', difftime(endTime[3], startTime[3], units = 'mins'), 'mins')
  cat('\n\n')
}
## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'auto'))
######### 


######### Save the results to file ----------
## save training regions
writeVector(regions, file.path('output', fp_base, 'modellingData/ModelRegion/yearly_modRegions.gpkg'), overwrite = T)

## save the training data
saveRDS(modData, file = file.path('output', fp_base, 'modellingData/yearly_modData.Rdata'))

## save the rasters
writeRaster(allModels, filename = file.path('output', fp_base, 'models/rast/yearly_predictions.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)

## save the variable importance
write.csv(allVarimp, file.path('output', fp_base, 'models/evaluation/yearly_varimp.csv'), row.names = F)

## save response curves
write.csv(allRc, file.path('output', fp_base, 'models/evaluation/yearly_response_curves.csv'), row.names = F)

## save models averages
writeRaster(allFavAVG, filename = file.path('output', fp_base, 'models/rast/yearly_mods_avg_fav.tif'), 
            gdal = c('COMPRESS=DEFLATE'), overwrite = T)
######### 

## offload memory
rm(list = ls())
gc()
