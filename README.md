# TimeVaryingCoefficientJM
## Time-varying coefficient joint model of longitudinal and survival data

This repository includes the code for a Bayesian joint model that allows a time-varying coefficient to link the longitudinal and the survival processes, using P-splines.

Specifically:
* "**Simulate**": includes the simulation of a dataset.
* "**ModelJAGS**": includes the joint model for jags.
* "**PrepareData**": includes the preparation of the data.
* "**Fit**": includes the main code. Specifically, it loads the packages, runs the simulated data, runs the code to prepare the data, creates the model, runs the model in JAGS and saves the results.

How does it work:
* Download all files and place them in one folder.
* Set as working directory in R this folder.
* Run the code in "Fit" for fitting the joint model.
