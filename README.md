# GitWDRORestrictedIntervention
This Git contains the R code and plots of the Master thesis "Estimation under Restricted Interventions via Wasserstein Distributionally Robust Optimisation". The methods are in the "src" folder. For instance, the Martingale WDRO estimator can be found in "src" folder.

To perform the outlier experiment with PH-Anchor-RWP run:
```
Rscript ./simulations/simulationAnchor.R
```
To run the simulation under latent confounding in Section 6.2 run:
```
Rscript ./simulations/simulationWDRO.R
```
Figures are located in the "results" folder within the "simulations folder. For instance, the regularisation path of the Martingale WDRO:
<p align="center">  
    <img src="./simulations/results/regupath.svg" width=75% height=75%>
</p> 
