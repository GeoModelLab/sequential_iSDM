#########################################################################################
############ Algorithms selection (a): Create the model region and apply VIF ############
#########################################################################################

library(terra)
library(fuzzySim)
library(tidyverse)

## load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
world <- geodata::world(path = 'utils', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica

occs <- read.csv(file.path('output', fp_base, 'modellingData/all_occs.csv'))
occs <- cleanCoords(occs, coord.cols = c('long', 'lat'), uncert.col = NULL, abs.col = NULL)
occs <- vect(occs, geom = c('long', 'lat'), crs = crs(world))

plot(world, lwd = 0, main = 'GBIF presence occurrences')
plot(occs, col = 'red', cex = 0.3, add = T)
plot(world, lwd = 0.5, add = T)

predictors <- rast(file.path('output', fp_base, 'modellingData/predictors.tif'))
plot(predictors[[1:ceiling(nlyr(predictors)/2)]], legend = F)
plot(predictors[[(ceiling(nlyr(predictors)/2)+1):nlyr(predictors)]], legend = F)
######### 


######### Modelling region ----------
## create a buffer of diameter wd around each point to be used as study region
w <- width(aggregate(occs))   # compute the width of the occurrences
scaling <- 0.05               # scaling factor for the width
wd <- w * scaling             # diameter around each point
## save the mult object for later use
saveRDS(wd, file.path('output', fp_base, 'modellingData/ModelRegion/region_mult.rds'))

## compute the buffer
b <- buffer(occs, width = wd)
## get the region
region <- aggregate(b)
region <- crop(region, world)
plot(region, col = 'gold', main = 'Model training region', border = NA)
plot(world, add = T)
points(occs, cex = 0.4)
## save to disk
writeVector(region, file.path('output', fp_base, 'modellingData/ModelRegion/CV_modRegion.shp'), overwrite = T)

## crop predictors according to model region
predictorsCropped <- crop(predictors, region, mask = T, touches = F)
plot(predictorsCropped[[1:ceiling(nlyr(predictorsCropped)/2)]], legend = F)
plot(predictorsCropped[[(ceiling(nlyr(predictorsCropped)/2)+1):nlyr(predictorsCropped)]], legend = F)
## save to disk
writeRaster(predictorsCropped, filename = file.path('output', fp_base, 'modellingData/crossValidation/CV_predictorsCropped.tif'),
            filetype = 'GTiff', datatype = 'FLT4S', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
gc()
######### 


######### Data Gridding and Absence Selection ----------
## presence VS. absence of record
gridData <- gridRecords(rst = predictorsCropped, pres.coords = occs)
nrow(gridData)
## clean from NAs in predictors if any
nas <- gridData |> 
  filter(if_any(everything(), is.na))
if (nrow(nas) > 0) {
  nasSPAT <- vect(nas, geom = c('x', 'y'), crs = crs(world))
  gridData <- drop_na(gridData)
  ## check which is NA with an interactive plot
  plet(predictorsCropped[['dem']]) |>
    points(occs, col = 'red') |>
    points(nasSPAT)
}
table(gridData$presence)
## save the raw gridded data to csv
write.csv(data.frame(gridData), file.path('output', fp_base, 'modellingData/crossValidation/CV_griddedData.csv'), row.names = F)

## download a bias layer to be used in pseudo-absence selection
bias <- geodata::footprint(year = 2009, path = 'utils')
bias <- project(bias, predictors)
bias <- crop(bias, region, mask = T)
plot(bias)

## select (biased) pseudo-absences from the pool of 0s in the study area
mult_p = 10 # means that the presence-to-absence ratio will be of 1:10
## construct the modelling data
modData <- selectAbsences(gridData,
                       sp.cols = "presence", coord.cols = c("x", "y"), CRS = crs(occs),
                       mult.p = mult_p, seed = 1312, df = T, bias = bias)
table(modData$presence)
## save to disk
write.csv(data.frame(modData), file.path('output', fp_base, 'modellingData/crossValidation/CV_modData.csv'), row.names = F)
######### 


######## Predictors analysis (VIF + pairwse correlation) --------
vars <- collinear::collinear(modData, responses = 'presence', predictors = names(modData)[-c(1:4)],
                             max_cor = 0.75, max_vif = 5)
write.csv(data.frame(predictors = vars$presence$selection), file.path('output', fp_base, 'modellingData/selectedVariables.csv'))
vars$presence$selection
######### 

## offload memory
rm(list = ls())
gc()

