# Multi-period forecasting

This repository contains the implementation of the experiments and methods for *"Smooth multi-period forecasting with application to prediction of COVID-19 cases"* paper by Tuzhilina et al.

## Data

The data used for the COVID-19 experiments was downloaded from the COVIDcast API https://cmu-delphi.github.io/delphi-epidata/api/covidcast.html. 
* train_non-missing.rds - train dataset without unobserved values
* train_missing.rds - train dataset with unobserved values
* test.rds - test dataset


## Code

The code includes three files.
* functions.R - the file containing the implementation of baseline MPF, smooth MPF, and QMPF, as well as several functions for graphical visualization of the results
* simulation.R - the code containing the simulation experiments described in the paper
* realdata.R - the implementation of the real data CovidCast experiments
