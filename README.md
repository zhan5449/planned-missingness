This is the accompanying Github repository to the paper 
"Planned Missngness to Reduce Survey Length: A Sheep in Wolf's Clothing" 
by Charlene Zhang, Paul Sackett, and Saron Demeke. 
We conducted a Monte Carlo simualtion that compares a planned missingness design 
(with multiple imputation) with using short forms for the goal of reducing survey length. 

This repository contains our scripts for conducting the simulation and 
analyzing the simulation results. You can also find the full simulated data, 
as well as supplementary analyses in addition to the results reported in the paper.

# Set up
1. Run Dockerfile to install needed R packeges.
(Windows: choco install -y docker-desktop; OS X: brew install --cask docker; 
Linux: Follow steps described in: <a href="https://docs.docker.com/engine/install/linux-postinstall/">Post-installation steps for Linux</a></td>)
2. Open an R console or RStudio window.
(R can be downloaded for free from https://cran.r-project.org; RStudio can be downloaded for free from https://www.rstudio.com/)

# Simulation
1. Run <em>population generation.R</em> to generate <em>Population Correlation Matrices.rds</em>, 
which contains the population correlation matrices on which simulations are based.
2. Run <em>SF vs PM simulation.R</em> to run entire simulation study. As the actual simulation 
was conducted via supercomputing resources, we provide the entirety of the simulation output, 
which can be found in the respective subfolders in the <em>data/ </em>folder.
3. Run <em>SF vs PM analysis.R</em> to process the simulation output and conduct the analyses 
included in the paper. Results are saved in the <em>analysis_results/</em> folder. We also tested
regression assumption for the regression models that we ran. The diagnostic plots can be found in the 
<em>distribution_diagnostics/</em> folder.
4. [supplementary materials] In an effort to understand the impact of the MICE settings that we used
for model convergence, we conducted detailed diagnostics for select conditions for manual inspection.
<em>simulation2_diagnostics.R</em> conducts the diagnostics and saves the output in the <em>two-factor diagnostics/</em> folder. 
<em>simulatino2_diagnostics_plots.R</em> takes the output and creates diagnostics plots. We also compared the distributions of the imputed values with 
that of the defined populations using <em>SF vs PM diagnostics_imputedvobserved.R</em>, which saves 
density plots in the <em>density_plots/</em> folder. 


