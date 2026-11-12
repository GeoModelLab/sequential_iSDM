# Sequential iSDM
This repository stores the code used in the publication: **Integrating the temporal dimension in invasive species distribution modelling: a case study on *Halyomorpha halys*** (DOI: <tbc>) that introduces the sequential iSDM pipeline to produce species distribution models for invasive organisms.

## Citation
Rossi et al. (in prep.) Integrating the temporal dimension in invasive species distribution modelling: a case study on Halyomorpha halys [Journal TBC]

## Overview
#### Article Abstract:
> Species Distribution Models of invasive species (iSDMs) are widely used to identify potentially suitable areas for alien species establishment. Standardly, they are developed using the entire set of available observations of the target organism, without accounting for the temporal evolution of its geographical spread. Here, we developed a sequential iSDM that iteratively evaluates the species environmental niche during its invasion process, using Halyomorphahalys, i.e., the Brown Marmorated Stink Bug (BMSB), as a case study. Results show that the performance of the sequential iSDMs largely changed along the simulation timeline, with minimum accuracy between 2008 and 2011 and a progressive convergence to the performance of traditional iSDMs onwards. This evidence supports the hypothesis of changes in the species realized ecological niche throughout the invasion process. Niche analyses explained these shifts and identified four distinct phases in BMSB invasion: evasion (2002 to 2008), expansion (2009 to 2015), consolidation (2016 to 2018), and saturation (2019 to 2025). Iterative model training addresses the violation of the equilibrium assumption in iSDMs, and our results show that the temporal dimension needs to be built explicitly into these models. This perspective can improve our understanding of past invasion trajectories and future establishment potential, supporting more informed risk reduction strategies in cultivated and natural ecosystems.

<br>

The sequential iSDM is defined as follows: let $\mathcal{D}_N$ be the set of species occurrences within the native range; $\mathcal{D}_I$ the set of occurrences from the invaded range; and $\mathcal{D}_I^{(t)}$ the invasive occurrences recorded until year $t$. For a given year $y$, the training dataset of the sequential iSDM is

$$
D_y = \mathcal{D}_N \cup \bigcup_{t=t_0}^y \mathcal{D}_I^{(t)} 
$$

where $t_0$ is first year of the invasion time series. Then, for each year $y$ a sequential iSDM is trained as 

$$
f_y(x) = \mathcal{M}(D_y, x)
$$

where $x$ is the set of predictors and $\mathcal{M}$ represents the iSDM fitting procedure. The resulting sequence

$$
\{f_{t_0}, f_{t_0+1}, ..., f_T\}
$$

provides a temporal representation of how the estimated environmental niche changes as new regions are invaded.

## Installation
### 1. Clone the repository
To clone this repository run the following code in a terminal in the desired location.
```
git clone git@github.com:camalloo/sequential_iSDM.git
```
### 2. Prepare the environment
Source the `00_Environment.R` script to create all required folders and download the necessary libraries.

## Required Data
GBIF occurrence data can be downloaded at the following [link](https://www.gbif.org/occurrence/search?taxonKey=PJGX6). The paper uses data up to 2025 for training and 2026 data for validation purposes only. Inside `data/observations` are stored two folders, one for training and one for validation data, where downloaded observations should be stored.  
Predictor variables should be stored in `data/predictors` as .tif files. The script contains a pipeline to automatically download the variables selected in the study from [CHELSA](https://www.chelsa-climate.org/) and saves them in `data/predictors/raw`.

## Code
The repository contains eight `.R` files inside the `scripts/` folder to reproduce the results of the paper. The code is entirely built in R (version 4.6.0):
- `00_Environment.R` sets up the correct working environment by creating the required folder structure and downloading the necessary R packages
- `01_Observations.R` processes the raw observation file downloaded from GBIF
- `02_Predictors.R` downloads and processes the predictor variables
- `03a` to `03c_Algorithms_selection.R` test SDM algorithms through a block cross validation procedure to select the best performing ones with respect the the data used
- `04_Static_models.R` computes the static iSDMs (all observations, native observations, invasive observations) using the best performing algorithms of 03
- `05_Yearly_runs.R` computes the sequential iSDMs using the best performing algorithms of 03
- `06_Evaluation.R` evaluates both static and sequential iSDMs using the Continuous Boyce Index and ROC-AUC metrics
- `07_Mess_PCA_Ensemble.R` performs the Multivariate Environmental Similarity Surface analysis, the PCA, and computes the Ensemble of models predictions

## Directory Tree
```
.
├── data
│   ├── observations
│   │   ├── training
│   │   └── validation
│   └── predictors
│       ├── raw
│       └── bioclim_chelsa.txt
├── scripts
│   └── ...
├── utils
│   └── assets
│       └── ...
├── LICENSE.txt
└── README.md
```

## Funding
The present project is supported by the Luxembourg National Research Fund (FNR), grant number: 19335551.

## License
Attribution-NonCommercial-NoDerivatives 4.0 International

<br>

<p align="center">
  <img src="utils/assets/logo.png" width="1000" alt="Logos" />
</p>
