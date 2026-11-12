#####################################################
############ Explore and analyse results ############
#####################################################

library(tidyverse)
library(terra)
library(fuzzySim)
library(pROC)
library(modEvA)
library(patchwork)

## load the mother folder
fp_base <- readRDS(file.path('utils/motherFolderPath.rds'))

######### Data ----------
world <- geodata::world(path = 'utils/', resolution = 1)
world <- world[!world$GID_0 == 'ATA'] # world map no Antactica

allObs <- read.csv(file.path(file.path('output', fp_base, 'modellingData/all_occs.csv')))
native <- read.csv(file.path('output', fp_base, 'modellingData/native_occs.csv')) |> 
  mutate(year = year(eventDate), status = 'n')
invasive <- read.csv(file.path('output', fp_base, 'modellingData/invasive_occs.csv')) |> 
  mutate(year = year(eventDate), status = 'i') |> filter(year >= 1996) # 1996 is first invasive discovery

## validation data are from the most recent year (not included in any training)
validationRaw <- read.csv('data/observations/validation/validation_occurrences.csv')

## all observations models
avgAllObs <- rast(file.path('output', fp_base, 'models/rast/AO_mods_avg_fav.tif'))
allObsVarimp <- read.csv(file.path('output', fp_base, 'models/evaluation/AO_varimp.csv'))
allObsRc <- read.csv(file.path('output', fp_base, 'models/evaluation/AO_response_curves.csv'))

## native models
avgNative <- rast(file.path('output', fp_base, 'models/rast/NT_mods_avg_fav.tif'))
nativeVarimp <- read.csv(file.path('output', fp_base, 'models/evaluation/NT_varimp.csv'))
nativeRc <- read.csv(file.path('output', fp_base, 'models/evaluation/NT_response_curves.csv'))

## invasive models
avgInvasive <- rast(file.path('output', fp_base, 'models/rast/IV_mods_avg_fav.tif'))
invasiveVarimp <- read.csv(file.path('output', fp_base, 'models/evaluation/IV_varimp.csv'))
invasiveRc <- read.csv(file.path('output', fp_base, 'models/evaluation/IV_response_curves.csv'))

## yearly models
avgYearly <- rast(file.path('output', fp_base, 'models/rast/yearly_mods_avg_fav.tif'))
yearlyVarimp <- read.csv(file.path('output', fp_base, 'models/evaluation/yearly_varimp.csv'))
yearlyRc <- read.csv(file.path('output', fp_base, 'models/evaluation/yearly_response_curves.csv'))
yearlyVarimp$Variable <- gsub('hursmean', 'hurs', yearlyVarimp$Variable)
yearlyRc$var <- gsub('hursmean', 'hurs', yearlyRc$var)
######### 


######### Create Validation Data ----------
## clean independent validation presences
valClean <- validationRaw |> 
  dplyr::select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, countryCode,
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
valClean <- valClean |> 
  CoordinateCleaner::cc_dupl(lon = 'decimalLongitude', lat = 'decimalLatitude', addition = 'eventDate', value = 'clean') |> 
  CoordinateCleaner::clean_coordinates(lon = 'decimalLongitude', lat = 'decimalLatitude', countries = 'countryCode', 
                                       tests = c('centroids', 'equal', 'gbif', 'outliers', 'seas', 'zeros'),
                                       value = 'clean') |> 
  dplyr::select(gbifID, occurrenceID, basisOfRecord, occurrenceStatus, species, continent, eventDate, year, month, day, countryCode,
         lat = decimalLatitude, long = decimalLongitude)
## additional cleaning with FuzzySim
valClean <- fuzzySim::cleanCoords(valClean, coord.cols = c('long', 'lat'), uncert.col = NULL)
valClean <- valClean |> 
  dplyr::select(eventDate, lat, long)

## sample independent background points where no presence was recorded
spatValidation <- vect(valClean, geom = c('long', 'lat'), crs = crs(world))
totalPres <- rbind(valClean, allObs)
totalPres <- vect(totalPres, geom = c('long', 'lat'), crs = crs(world))

## absence points on the world
stackedEns <- c(avgAllObs[['mods_mean_fav_']], avgYearly[['mods_mean_fav_2025']],
                avgNative[['mods_mean_fav_']], avgInvasive[['mods_mean_fav_']]) ## stack the ensembles
bg <- gridRecords(rst = stackedEns, pres.coords = totalPres)
bg <- bg |>
  filter(presence == 0) |> 
  dplyr::select(presence, x, y)
table(bg$presence)
bgSPAT <- vect(bg, geom = c('x', 'y'), crs = crs(world))
## plot
plot(world, main = 'Validation Dataset', font.main = 1)
points(bgSPAT, col = 'blue', cex = 0.5)
points(spatValidation, col = 'red', cex = 0.5)
#########  


######### Validate Yearly Models ----------
## prepare observation dataframes for validation
valInvasive <- native |>  # invasive are validated with native
  mutate(presence = 1) |> 
  dplyr::select(presence, x = long, y = lat)
valNative <- invasive |> # native are validated with invasive
  mutate(presence = 1) |> 
  dplyr::select(presence, x = long, y = lat)
valInd <- valClean |>
  mutate(presence = 1) |> # most recent year data used in all validations
  dplyr::select(presence, x = long, y = lat)

## run validation script
iter <- 10                            # define the number of iteration for the validation
allEvals <- data.frame()              # store the results
# sampleSize <- min(nrow(valInvasive),  # use the same number of presences for all models in validation
#                   nrow(valNative), 
#                   nrow(valInd))
## go through all years
years <- names(readRDS(file.path('output', fp_base, 'modellingData/yearly_modData.Rdata')))
set.seed(1312) # seed for random sampling of validation presences and absences
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  for (yr in 1:length(years)) {
    y <- years[yr]
    cat('\t\tEvaluating model year', y, '\n')
    
    ## get the validation set
    if (y == max(years)) { # if last year of training evaluate with most recent records not used fro training (valClean)
      val <- valInd
    } else {
      val <- invasive |> 
        filter(year > y) |> 
        mutate(presence = 1) |> 
        dplyr::select(presence, x = long, y = lat)
      val <- rbind(val, valInd)
    }
    valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
    
    ## 1. handle absences maintaining a 1:1 prevalence ratio
    # bgrnd <- bg |> slice_sample(n = sampleSize)       # OR
    ## 2. handle absences maintaining a 10:1 prevalence ratio
    bgrnd <- bg |> slice_sample(n = nrow(val)*10)
    bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
    
    ## get the ensemble prediction for this year to be evaluated
    toEval <- avgYearly[[paste0('mods_mean_fav_', y)]]
    
    ## extract prediction at validation set
    prediction <- gridRecords(rst = toEval, pres.coords = valSPAT, abs.coords = bgrndSPAT)
    names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
    
    ## AUC
    roc <- roc(response = prediction$presence, 
               predictor = prediction$pred, 
               quiet = TRUE)
    auc <- as.numeric(auc(roc))
    auc_ci <- as.numeric(ci.auc(roc))  # lower, mean, upper
    
    ## Boyce index
    boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
    
    ## store results
    allEvals <- rbind(allEvals,
                      data.frame(iteration = idx,
                                 tag = 'yearly',
                                 year = y,
                                 auc = auc,
                                 auc_lower = auc_ci[1],
                                 auc_upper = auc_ci[3],
                                 boyce = boyce))
  }
}
#########  


######### Validate AO Model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- valInd
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgAllObs[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  
  ## store results
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'all_obs_model',
                               year = 2,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce))
}
#########


######### Validate NT Model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- rbind(valNative, valInd)
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgNative[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  
  ## store result
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'native_model',
                               year = 0,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce))
}
#########


######### Validate IV Model ----------
set.seed(1312)
for (idx in 1:iter) {
  cat('iteration', idx, '\n')
  
  ## validation presences
  val <- rbind(valInvasive, valInd)
  valSPAT <- vect(val, geom = c('x', 'y'), crs = crs(world))
  
  ## handle absences maintaining a 10:1 prevalence ratio
  bgrnd <- bg |> slice_sample(n = nrow(val)*10)
  bgrndSPAT <- vect(bgrnd, geom = c('x', 'y'), crs = crs(world))
  
  ## extract prediction at validation presences
  prediction <- gridRecords(rst = avgInvasive[['mods_mean_fav']], pres.coords = valSPAT, abs.coords = bgrndSPAT)
  names(prediction) <- c('presence', 'x', 'y', 'cells', 'pred')
  
  ## AUC
  roc_all <- roc(response = prediction$presence, 
                 predictor = prediction$pred, 
                 quiet = TRUE)
  auc_all <- as.numeric(auc(roc_all))
  auc_all_ci <- as.numeric(ci.auc(roc_all))  # lower, mean, upper
  
  ## Boyce
  boyce <- Boyce(obs = prediction$presence, pred = prediction$pred, plot = F)$B
  
  ## store results
  allEvals <- rbind(allEvals,
                    data.frame(iteration = idx,
                               tag = 'invasive_model',
                               year = 1,
                               auc = auc_all,
                               auc_lower = auc_all_ci[1],
                               auc_upper = auc_all_ci[3],
                               boyce = boyce))
}

## aggregate all evals
allEvalsAgg <- allEvals |> 
  group_by(tag, year) |> 
  summarise(auc_mean = mean(auc),
            auc_sd = sd(auc),
            boyce_mean = mean(boyce),
            boyce_sd = sd(boyce))
write.csv(allEvalsAgg, 
          file.path('output', fp_base, 'models/evaluation/final_evaluation_all_models.csv'),
          row.names = F)
#########  


######### Plot performance metrics ----------
allEvalsAgg <- read.csv(file.path('output', fp_base, 'models/evaluation/final_evaluation_all_models.csv'))
years <- unique(allEvalsAgg$year)

## CBI
cbiPLOT <- ggplot(allEvalsAgg |> dplyr::filter(tag == 'yearly' & !year %in% c(0, 1, 2)),
                  aes(x = year)) +
  geom_ribbon(aes(y = boyce_mean, ymin = boyce_mean-boyce_sd, ymax = boyce_mean+boyce_sd), 
              colour = 'blue', fill = 'blue', linewidth = 0.1, alpha = 0.4) +
  geom_line(aes(y = boyce_mean), colour = 'blue') +
  # geom_hline(aes(yintercept = 1), linetype = 'solid', colour = 'black') +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 2], linetype = 'AO')) +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 0], linetype = 'NT')) +
  geom_hline(aes(yintercept = allEvalsAgg$boyce_mean[allEvalsAgg$year == 1], linetype = 'IV')) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_linetype_manual(values = c('AO' = 'dotdash',
                                   'NT' = 'dashed',
                                   'IV' = 'dotted'),
                        breaks = c('NT', 'IV', 'AO'),
                        labels = c('AO' = 'All Obs. (AO)',
                                   'NT' = 'Native (NT)',
                                   'IV' = 'Invasive (IV)')) +
  labs(x = 'Training Year',
       y = 'Metric Score',
       title = '(A) CBI',
       linetype = 'Static SDMs') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

## AUC
aucPLOT <- ggplot(allEvalsAgg |> dplyr::filter(tag == 'yearly' & !year %in% c(0, 1, 2)),
                  aes(x = year)) +
  geom_ribbon(aes(y = auc_mean, ymin = auc_mean-auc_sd, ymax = auc_mean+auc_sd), 
              colour = 'red', fill = 'red', linewidth = 0.1, alpha = 0.4) +
  geom_line(aes(y = auc_mean), colour = 'red') +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 2], linetype = 'AO')) +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 0], linetype = 'NT')) +
  geom_hline(aes(yintercept = allEvalsAgg$auc_mean[allEvalsAgg$year == 1], linetype = 'IV')) +
  scale_x_continuous(breaks = years) +
  scale_y_continuous(limits = c(0.3, 1)) +
  scale_linetype_manual(values = c('AO' = 'dotdash',
                                   'NT' = 'dashed',
                                   'IV' = 'dotted'),
                        breaks = c('NT', 'IV', 'AO'),
                        labels = c('AO' = 'All Obs. (AO)',
                                   'NT' = 'Native (NT)',
                                   'IV' = 'Invasive (IV)')) +
  labs(x = 'Training Year',
       y = 'Metric Score',
       title = '(B) AUC',
       linetype = 'Static SDMs') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1))

evaluationPLOTS <- cbiPLOT+aucPLOT + plot_layout(guides = 'collect', axes = 'collect')
evaluationPLOTS
######### 


######### Visualize Varimp and Response Curves ----------
## combine response curves and variable importance of all models in a single plot
## response curves
ggplot() +
  geom_line(data = yearlyRc, aes(x = x, y = y, color = factor(year))) +
  geom_line(data = baselineRc, aes(x = x, y = y, linetype = year), color = 'black') +
  scale_linetype_manual(
    name = 'Static SDMs',
    values = c( 'AO' = 'solid', 'IV' = 'dotted', 'NT' = 'dashed')) +
  facet_grid(mod~var, scales = 'free') +
  labs(title = '(A) Partial Dependence Plots',
       x = 'Variable',
       y = 'Prediction',
       color = 'Dynamic SDMs') +
  guides(colour = guide_legend(ncol = 1)) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## variable importance
years <- sort(unique(na.omit(yearlyVarimp$year))) |> as.character()
combinedVarimp <- bind_rows(baselineVarimp |> 
                              mutate(progression = as.factor(tag)), 
                            yearlyVarimp |> 
                              mutate(progression = as.factor(year)))
combinedVarimp <- combinedVarimp |> 
  mutate(progression = as.character(progression)) |> 
  mutate(progression = factor(progression,
                              levels = c('NT', years, 'IV', 'AO'))) |>
  arrange(progression)

ggplot(data = combinedVarimp, aes(x = progression, y = Importance)) +
  geom_line(aes(group = Variable, color = Variable)) +
  facet_wrap(~mod, nrow = 3) +
  labs(title = 'Variable Importance Plot',
       x = 'Year',
       y = 'Importance',
       color = 'Variable') +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
######### 

## offload memory
rm(list = ls())
gc()
