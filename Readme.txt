### README b_Bayesian_2024
## Probabilistic estimate of b-value temporal variations
## Update: 13/11/24
## Author: Marine Laporte
## Contact: marinelaporte@icm.csic.es
## Citation: 
###-----------------------------------------

This Readme describes how to use b-Bayesian to perform the transdimensional inversion of the temporal evolution of :
	- the b-value (from the Gutenberg-Richter law)
	- mu and sigma (describing the detectability as defined by the error function from Ogata and Katsura (1993))	
	
The code returns : 
	- the probability of the number and location of temporal discontinuities : prop.txt (location) / Nc_tot (number)
	- the full probability distribution of b-value over time : Likelihood_b.txt
	- the full probability distribution of mu over time : Likelihood_mu.txt
	- the full probability distribution of sigma over time : Likelihood_sigma.txt

> ⚠️ **This is not the latest version of b-Bayesian.**
> This is the original version described in Laporte et al. (2025), kept for reproducibility of the published results.
> For any new study, use the latest and recommended version:
> 👉 [`b-Bayesian_binned_McMC_2026`](../b-Bayesian_binned_McMC_2026) – full McMC, faster, accounts for magnitude discretization.

#############################
##  INSTALLATION
#############################

# Dependancies 
	- for running b-Bayesian : fortran90 compiler (to choose in the Makefile)
	- for creating figures with the outputs : python3 
		libraries : numpy, matplotlib, crameri, math
		
#############################
##  Inversion b-Bayesian
#############################

# Change input file/ parametrisation :
In the script /inversion/b_Bayesian.f90

	- Change data path : Magnitude_file = "your data path"
	- Change parameters of the Markow chains (number of iterations to reach convergence, burn-in ...)
	- Adapt the priors for beta, mu, sigma to your dataset (large priors to begin with)


# Run b-Bayesian
In the folder /inversion/ : > make	
		            > ./run
		            
#############################
##  OUTPUTS
#############################

b-Bayesian returns 8 output text files in the /outputs folder :
	- input.txt : a copy of the earthquake dataset
	- hypercube.txt : a save of the hyper-parameters chosen for the run (priors, number of iterations) 
	- prop.txt (size : T_bins) : the final probability of a temporal discontinuity for for each temporal bins of T_bins
	- Nc_tot.txt (size : itmax) : the number of temporal discontinuities at each iteration of the McMc -> which converges theoretically towards the target number 
	- Likelihood_b (size : T_bins x B_bins) : Marginal probability density of b-value over time 
	- Likelihood_mu (size : T_bins x Mu_bins) : Marginal probability density of mu over time 
	- Likelihood_sigma (size : T_bins x Sigma_bins) : Marginal probability density of sigma over time 
	- Acceptances.txt (3x3): Columns = [Birth Death Moves] // Lines = [Number of accepted models ; Number of proposed models ; Final acceptance rate of the births/deaths/moves ]


#############################
##  PLOTS FIGURES
#############################

In /figures the python code output_figures.py reads the outputs generated in the folder /outputs and creates several figures :


# Main :

- Rename your output-figures in the main of output_figures.py

> python3 output_figures_bMuSigma.py

# Returns 

Creates a new folder in the /figures folder : /Run_day_month_year_hour

Fig0 : Plots the magnitudes over times and surimposed the probability of discontinuities
Fig1 : 4 subplots : (up-left) Final number of accepted temporal discontinuities for each T_bins
		    (up-right) Number of accepted models vs Number of proposed models for each cases birth/death/move
		    (down-left) Number of discontinuities/model at each iteration
		    (down-right) Histogram associated to the down-left figure -> The final number of discontinuities is the max of the histogram
Fig2 : Marginal probability density over time for b-value
Fig3 : Marginal probability density over time for mu
Fig4 : Marginal probability density over time for sigma


