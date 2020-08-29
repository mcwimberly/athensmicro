# athensmicro

This repository contains the code used to carry out the analyses and generate the graphics for this paper.

Wimberly, M. C., J. K. Davis, M. V. Evans, A. Hess, P. M. Newberry, N. Solano-Asamoah, C. C. Murdock. In Press. Land cover affects microclimate and temperature suitability for arbovirus transmission in an urban landscape. PLoS Neglected Tropial Diseases.

proc_microdata.R: Initial processing of microclimate data to calculate minimum and maximum daily temperatures.

downscale_temp.R: Generate 30 m resolution microclimate maps by combining microclimate data with land cover maps and macroclimate grids.

fit_abundance.R: Fit an empirical model of *Aedes albopictus* abundance as a function of minimum and maximum temperature.

calculate_vc.R: Estimate the vectorial capacity using a mechanistic model based on temperature dependent mosquito traits.

generate_manuscript_figures.R: Contains the code used to generate all figures in the manuscript.


