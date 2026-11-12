######################################################
############ Process the raw observations ############
######################################################

library(tidyverse)
library(terra)

## load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
## download world map and store it for future use in utils
world <- geodata::world(path = 'utils', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antarctica

## get the observation data
data_fp <- 'data/observations/training/occurrences.csv'
dataRaw <- read.csv(data_fp)
######### 


######### Cleaning ----------
## clean the data
dataClean <- dataRaw |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
         decimalLatitude, decimalLongitude, verbatimEventDate, year, month, day, startDayOfYear, endDayOfYear) |> 
  mutate(eventDate = as.Date(eventDate, '%d/%m/%Y')) |> 
  ## make the date from day, month, year when available
  mutate(eventDate = coalesce(eventDate, as.Date(paste(day, month, year, sep = '/'), '%d/%m/%Y'))) |> 
  ## make the date from doy and year when available 
  mutate(eventDate = coalesce(eventDate, as.Date(startDayOfYear-1, origin = paste0(year, '/01/01')))) |>
  ## convert verbatim dates
  mutate(eventDate = coalesce(eventDate, as.Date(verbatimEventDate, format = '%d/%m/%Y %H:%M'))) |> 
  ## drop no coords or date
  drop_na(decimalLatitude, decimalLongitude, eventDate)

## remove duplicated records and check for other error in coordinates --> FINAL DF
dataClean <- dataClean |> 
  CoordinateCleaner::cc_dupl(lon = 'decimalLongitude', lat = 'decimalLatitude', addition = 'eventDate', value = 'clean') |> 
  CoordinateCleaner::clean_coordinates(lon = 'decimalLongitude', lat = 'decimalLatitude', countries = 'countryCode', 
                                       tests = c('centroids', 'equal', 'gbif', 'outliers', 'seas', 'zeros'),
                                       value = 'clean') |> 
  select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, year, month, day, countryCode,
         lat = decimalLatitude, long = decimalLongitude)
## additional cleaning with FuzzySim
dataClean <- fuzzySim::cleanCoords(dataClean, coord.cols = c('long', 'lat'), uncert.col = NULL)

## plot all
plot(world)
points(dataClean$long, dataClean$lat, pch = 20, cex = 0.3, col = 'red')
######### 


######### Split Native from Invasive and Save ----------
## occurrences from the native areas
native <- dataClean |> 
  filter(countryCode %in% c('CN', 'TW', 'JP', 'KP', 'KR', 'HK', 'MO'))
## plot
plot(world, main = expression(paste(italic('H. halys'), ' ditribution')))
points(native$long, native$lat, col = 'blue', pch = 16, cex = 0.3)

## occurrences from invaded areas
invasive <- dataClean |> 
  filter(!(countryCode %in% c('CN', 'TW', 'JP', 'KP', 'KR', 'HK', 'MO')))
## add to plot
points(invasive$long, invasive$lat, col = 'red', pch = 16, cex = 0.3)

## legend of the distribution plot
legend('bottomleft', horiz = F, legend = c('Native Range', 'Invaded Range'), 
       y.intersp = 0.2, xpd = T, inset = c(0, 0.15),
       col = c('blue', 'red'), pch = 16, pt.cex = 1, bty = 'n')

## save cleaned occs to .csv
write.csv(native |> select(eventDate, lat, long), file.path('output', fp_base, 'modellingData/native_occs.csv'), row.names = F)
write.csv(invasive |> select(eventDate, lat, long), file.path('output', fp_base, 'modellingData/invasive_occs.csv'), row.names = F)
write.csv(dataClean |> select(eventDate, lat, long), file.path('output', fp_base, 'modellingData/all_occs.csv'), row.names = F)
#########

## offload memory
rm(list = ls())
gc()
