# WDRORestrictedIntervention Estimation
This Git contains the R code and plots of the Master Thesis "Estimation under Restricted Interventions via Wasserstein Distributionally Robust Optimisation". The methods are in the "src" folder and the simulations are in the "simulations" folder. Within the "simulations" folder we have figures in the "results" folder.    

## Simulations
Section 6.1:

To perform the outlier experiment with PH-Anchor-RWP run:
```
Rscript ./simulations/simulationAnchor.R
```
To perform the anchor transformation experiment run:
```
Rscript ./simulations/simulationAnchorTransformation.R
```
Section 6.2:

To run the simulation under latent confounding in Section 6.2 run:
```
Rscript ./simulations/simulationWDRO.R
```
To obtain the plots run:
```
Rscript ./simulations/plot_generator.R
```
To obtain the regularisation path run:
```
Rscript ./simulations/regu_path.R
```
To do the Ancestor WDRO run:
```
Rscript ./simulations/simulationAncestorWDRO.R
```
Section 6.3:

To obtain the Wasserstein causal regulariser run:
```
Rscript ./simulations/simulationWassersteinCausalRegulariser.R
```

## Comments
Any feedback or comments on the code or methods is more than welcome. Correspondence is at [email](mailto:mvanden@student.ethz.ch).
