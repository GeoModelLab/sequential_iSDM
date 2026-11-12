##################################################################
############ Prepare the environment for the analysis ############
##################################################################


######### Install Packages ----------
install.packages("tidyverse", "terra", "fuzzySim",
                 "modEvA", "blockCV", "predicts",
                 "maxnet", "mgcv", "dbarts",
                 "randomForest", "embarcadero", "vip",
                 "pdp", "pROC", "patchwork",
                 "tidyterra", "FactoMineR", "factoextra",
                 "dismo", "pheatmap", )
source('https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df')
######### 


######### Create Folders ----------
## data/observations
dir.create('data/observations/training', recursive = T)
dir.create('data/observations/validation', recursive = T)
## data/predictors
dir.create('data/predictors/raw', recursive = T)

## output
dir.create('output')

## utils
dir.create('utils')

## output/fp_base --> where the result of a run are stored
fp_base <- 'Run_1'
## these subdirs must be included: fp_base/subdirs
subdirs <- c('ensemble',
             'mess',
             'modellingData',
             'modellingData/crossValidation',
             'modellingData/ModelRegion',
             'models',
             'models/crossValidation',
             'models/evaluation',
             'models/rast',
             'niche')
## the complete path 
path <- file.path('output', fp_base)
## create the subdirectories of output/fp_base
if (!dir.exists(path)) {
  dir.create(path)
  for (d in subdirs) {
    dir.create(file.path(path, d))
  }
  message('All output dirs created')
} else {
  ## check if all subdir exixt in the onther folder
  if (all(dir.exists(file.path(path, subdirs)))) {
    message('All output dirs exist already')
  } else {
    message('Missing some output dirs...')
    missingDirs <- file.path(path, subdirs)[!dir.exists(file.path(path, subdirs))]
    for (md in missingDirs) {
      dir.create(file.path(md), recursive = T)
    }
    message(paste('... all missing output dirs created'))
  }
  
}

## save the name of the mother folder as rds so stays consistent along the scripts
saveRDS(fp_base, file = file.path('utils/motherFolderPath.rds'))
######### 
