#############################################################################################
############ Run the models with all observation and project on the entire world ############
#############################################################################################

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

native <- read.csv(file.path('output', fp_base, 'modellingData/native_occs.csv'))
invasive <- read.csv(file.path('output', fp_base, 'modellingData/invasive_occs.csv'))
allObs <- read.csv(file.path(file.path('output', fp_base, 'modellingData/all_occs.csv')))

predictors <- rast(file.path('output', fp_base,'modellingData/predictors.tif'))
bias <- geodata::footprint(year = 2009, path = 'utils')
bias <- project(bias, predictors)

selectedMods <- read.csv(file.path('output', fp_base, 'models/crossValidation/selectedModels.csv'))
vars <- read.csv(file.path('output', fp_base, 'modellingData/selectedVariables.csv'))[, 2] # variables after collinearity analysis
## modelling region multiplier
wd <- readRDS(file.path('output', fp_base, 'modellingData/ModelRegion/region_mult.rds'))
######### 


######### Run Static Models ----------
## selected models
selectedMods$mod
## select the multiplier of absences
mult_p <- 10
## time tracking
(start <- Sys.time())
## loop over the three scenarios: native = NT; invasive = IV; all observations = AO
for (scenario in c('NT', 'IV', 'AO')) {
  cat('Scenario:', scenario, '\n')
  ## yearly time tracking
  startTime <- proc.time()
  
  ## modelling data ----
  if (scenario == 'NT') {
    occs <- vect(native, geom = c('long', 'lat'), crs = crs(world))
  } else if (scenario == 'IV') {
    occs <- vect(invasive, geom = c('long', 'lat'), crs = crs(world))
  } else if (scenario == 'AO') {
    occs <- vect(allObs, geom = c('long', 'lat'), crs = crs(world))
  }

  ## define the model region
  b <- buffer(occs, width = wd)
  region <- aggregate(b)
  region <- crop(region, world)
  writeVector(region, file.path('output', fp_base, 'modellingData/ModelRegion', paste0(scenario, '_modRegion.shp')), overwrite = T)
  ## crop predictors according to model region
  predictorsCropped <- crop(predictors, region, mask = T, touches = F)
  gc()
  
  ## data gridding
  gridData <- gridRecords(rst = predictorsCropped, pres.coords = occs)
  nrow(gridData)
  ## clean from NAs in predictors if any
  gridData <- drop_na(gridData)
  ## crop bias layer for absence generation to study area
  biasCropped <- crop(bias, region, mask = T)
  ## absences selection
  data <- selectAbsences(gridData,
                         sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                         mult.p = mult_p, seed = 1312, plot = F, df = T, bias = biasCropped, verbosity = 0)
  ## save modelling data to disk
  write.csv(data.frame(data), file.path('output', fp_base, 'modellingData', paste0(scenario, '_modData.csv')), row.names = F)
  #####
  
  ## model runs ----
  cat('\tGAM ')
  gamForm <- as.formula(paste0('presence ~ ', paste0('s(', vars, ', bs = "cs", k = 20)', collapse = '+')))
  gamMod <- gam(gamForm, family = binomial, data = data, model = T, method = 'REML')
  gamPred <- predict(predictors, gamMod, type = 'response', na.rm = T)
  names(gamPred) <- 'GAM.pred'
  gamFav <- Fav(pred = gamPred, sample.preval = prevalence(model = gamMod))
  names(gamFav) <- 'GAM.fav'
  gc()
  
  cat('- BART ')
  bartMod <- bart(x.train = data[, vars], y.train = data[, 'presence'],
                  ntree = 200, keeptrees = T, verbose = F, seed = 1312)
  bartPred <- predict2.bart(bartMod, raster::stack(predictors), quantiles = c(), splitby = 10, quiet = T)
  names(bartPred) <- c('BART.pred')
  bartPred <- rast(bartPred)
  bartFav <- Fav(pred = bartPred, sample.preval = prevalence(model = bartMod))
  names(bartFav) <- c('BART.fav')
  gc()
  
  cat('- RF\n')
  rfForm <- formula(paste0('presence ~ ', paste(vars, collapse = '+')))
  rfMod <- randomForest(rfForm, data, na.action = na.exclude)
  rfPred <- predict(predictors, rfMod, na.rm = T, cores = 10, cpkgs = 'randomForest')
  rfPred <- clamp(rfPred, lower = 0, upper = 1) # clamp to 0,1
  names(rfPred) <- 'RF.pred'
  rfFav <- Fav(pred = rfPred, sample.preval = prevalence(model = rfMod))
  names(rfFav) <- 'RF.fav'
  gc()

  ## stack predictions in a single raster
  all <- c(gamPred, gamFav,
           bartPred, bartFav,
           rfPred, rfFav)
  ## save to disk
  writeRaster(all, filename = file.path('output', fp_base, 'models/rast', paste0(scenario, '_predictions.tif')), 
              gdal = c('COMPRESS=DEFLATE'), overwrite = T)
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
  GAMvarimp <- cbind(mod = 'GAM', GAMvarimp)
  ## rf
  RFvarimp <- vi_permute(rfMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                         pred_wrapper = rf_pred_wrapper, nsim = 10,
                         parallel = T, .packages = 'randomForest')
  RFvarimp <- cbind(mod = 'RF', RFvarimp)
  ## bart
  BARTvarimp <- vi_permute(bartMod, train = subdata[, vars], target = subdata$presence, metric = 'RMSE', 
                           pred_wrapper = bart_pred_wrapper, nsim = 10,
                           parallel = T, .packages = 'dbarts')
  BARTvarimp <- cbind(mod = 'BART', BARTvarimp)
  
  ## save all to df
  ALLvarimp <- rbind(GAMvarimp, RFvarimp, BARTvarimp)
  write.csv(ALLvarimp, file.path('output', fp_base, 'models/evaluation', paste0(scenario, '_varimp.csv')))
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
               var = v,
               x = p[[v]],
               y = p$yhat)
  }) |> bind_rows()
  ## save all response curves to files
  ALLrc <- rbind(GAMrc, RFrc, BARTrc)
  write.csv(ALLrc, file.path('output', fp_base, 'models/evaluation', paste0(scenario, '_response_curves.csv')))
  cat('\tresponse curve done\n')
  #####
  
  ## models average ----
  favMean <- app(c(bartFav[[1]], gamFav, rfFav), mean)
  names(favMean) <- paste0('mods_mean_fav')
  favVar <- app(c(bartFav[[1]], gamFav, rfFav), var)
  names(favVar) <- paste0('mods_var_fav')
  favAVG <- c(favMean, favVar)
  writeRaster(favAVG, filename = file.path('output', fp_base, 'models/rast', paste0(scenario, '_mods_avg_fav.tif')), 
              gdal = c('COMPRESS=DEFLATE'), overwrite = T)
  #####
  
  ## scenario time tracking
  endTime <- proc.time()
  cat('Elapsed scenario', scenario, 'time:', difftime(endTime[3], startTime[3], units = 'mins'), 'mins')
  cat('\n\n')
}
## time tracking
(end <- Sys.time())
(length <- difftime(end, start, units = 'auto'))
######### 

## offload memory
rm(list = ls())
gc()
