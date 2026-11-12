########################################################
############ Download and select predictors ############
########################################################

library(tidyverse)
library(terra)

## load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
world <- geodata::world(path = 'utils', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica
######### 


######### Download BIOCLIM+ from CHELSA ----------
# links <- read.table('data/predictors/bioclim_chelsa.txt', header = FALSE)[,1]
# for (link in links) {
#   name <- basename(link)
#   download.file(link, destfile = file.path('data/predictors/raw', name))
# }
# 
# ## offload memory
# rm(link, links, name)
# gc()
######### 


######### Process Predictors ----------
## load vars
bio <- rast(list.files(paste0('data/predictors/raw'), full.names = T))
nlyr(bio) # check if all loaded
## handle names
names(bio) <- gsub('CHELSA_|_1981.*', '', names(bio))
names(bio)
## aggregate to 10km
bio <- aggregate(bio, fact = 10, fun = 'mean', na.rm = T, cores = 10)
res(bio)
## mask the world
bio <- crop(bio, world, snap = 'out', mask = T)

## plot them
plot(bio[[1:ceiling(nlyr(bio)/2)]], legend = F)
plot(bio[[(ceiling(nlyr(bio)/2)+1):nlyr(bio)]], legend = F)

# ## scale values to z-score (optional)
# stats <- global(bio, c('mean', 'sd'), na.rm = T)
# bio <- (bio-stats$mean) / stats$sd

## save to disk
writeRaster(bio, filename = file.path('output', fp_base, 'modellingData/predictors.tif'),
            filetype = 'GTiff', gdal = c('COMPRESS=DEFLATE'), overwrite = T)
######### 

## offload memory
rm(list = ls())
gc()
