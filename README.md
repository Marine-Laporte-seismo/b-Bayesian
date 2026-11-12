# Readme b-Bayesian
>author : Marine Laporte  
>contact : marine.laporte@univ-lyon1.fr  
>version : 28/08/2024  

> ⚠️ **This is not the latest version of b-Bayesian.**
> This is the original version described in Laporte et al. (2025), kept for reproducibility of the published results.
> For any new study, use the latest and recommended version:
> 👉 [`b-Bayesian_binned_McMC`](../b-Bayesian_binned_McMC_2026) – full McMC, faster, accounts for magnitude discretization.

This Readme describes how to use the Transdi_B code for performing a transdimensional inversion of the temporal evolution of  $\beta = \frac{b}{log(10)}$ from the Gutenberg Richter law and the two parameters $\mu$ and $\sigma$ of a given detection function law (default : error function). It is based on the fortran90 code RJ_MCMC of Thomas Bodin.  

This github combines f90 codes for performing the transdimensional inversion and two python codes for generating synthetic data and plotting the results of the transdimensional inversion.  

### 1 - Generate synthetic dataset /data
This branch contains a python code generate_synthetic_dist_etas.py. 

### 2 - Transdimensional inversion 
This branch contains the fortran 90 codes to proceed the transdimensional inversion. 
* main : Transdi_B.f90
* /src : likelihoods_random_marginals.f90
* /src : likelihoods_full_marginals.f90
* /out
### 3 - Plotting results after the inversion /visu
This branch contains a python code output_figures.py
The code reads all the outputs of /inversion/out 
* Figure 1 : Sum of temporal discontinuities
  Only a simple plot of the found discontinuities over time. Magnitudes over time are plotted in background to compare with the magnitudes distribution over time. 
  ![temporal_discontinuities](http://www.sampleurl.com/logo.png) 
* Figure 2 : Overview of the inversion
  A 4 subplots figure with :
  - up-left) same as figure 1
  - up-right) The % of acceptance for birth/death/pert
  - bot-left) The evolution of the number of discontinuities found with increasing iterations
  - bot-right) histogram of bot-left)
* Figure 3 : Likelihood mu
* Figure 4 : Likelihood b-value
* Figure 5 : Likelihood sigma
