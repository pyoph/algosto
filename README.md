## Overview


This repository contains the R code and supplementary materials required to reproduce the numerical experiments and real-data analyses presented in the accompanying manuscript.

The implementation relies on the STARRS R package, available from the Comprehensive R Archive Network (CRAN), together with other publicly available R packages.


## Repository organization

The supplementary materials are organized as follows:


├── README.md
├── main_SMD.R
├── main_simulated.R
├── R/
│   ├── algexec.R
│   ├── computeç_criteria.R
│   ├── crit_SMD.R
│   ├── parameters.R
│   └── plots.R
│   └── plots_SMD.R
│   └── simuls.R
│   └── smd_algo.R
│   └── smd_data_preprocessing
├── anonym_package/



## Required packages 

The required packages can be installed by running either main_SMD.R or main_simulated.R. In the anonymized package, an additional return value, time_cumulative, has been added to the onlineRobustVariance function.

## Experiments on the simulated dataset

The experiments on the simulated dataset can be run using main_simulated.R. The scripts parameters.R and simuls.R are used to generate the simulated dataset.

## Experiments on the real dataset

The real-data application uses the Server Machine Dataset (SMD).

The dataset is available from Kaggle: https://www.kaggle.com/datasets/mgusat/smd-onmiad . The data preprocessing steps can be performed by running the smd_data_preprocessing.R script.
 
