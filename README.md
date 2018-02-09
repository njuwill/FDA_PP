# Home Directory

This is a collection of codes and instruction documents to the project *Multilevel Functional Data Analysis for Temporal Point Processes with Applications in Stock Market Trading.*

## Project Summary

We propose a general framework of using multi-level log-Gaussian Cox processes to model repeatedly observed point processes with complex structures. A novel nonparametric approach is developed to consistently estimate the covariance kernels of the latent Gaussian processes at all levels. Consequently, multi-level functional principal component analysis can be conducted to investigate the various sources of variations in the observed point patterns. In particular, to predict the functional principal component scores, we propose a consistent estimation procedure by maximizing the conditional likelihoods of super-positions of point processes. We further extend our procedure to the bivariate point process case where potential correlations between the processes can be assessed. Asymptotic properties of the proposed estimators are investigated, and the effectiveness of our procedures is illustrated by a simulation study and an application to a stock trading dataset.

## Program

* Simulation: See `simulation` folder.
* Estimation   
   * Covariance Estimation and MFPCA. See `decomposition` folder.
   * Bandwidth Selection: h selection for the nonparametric covariance estimation. See `bandwidth selection` folder
* Prediction: Prediction of principal component scores. See `prediction` folder.
* MFPCA for Bivariate log-Gaussian Cox processes.   
  Here we extend the analysis of univariate case to bivariate case. See `bivariate estimation`.


