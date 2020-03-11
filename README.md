# ABA_functional_networks
The proc_all.m is a cell mode matlab file contains all steps
necessary to run analyses in sequential order. Set variables
in lines 5-8.The proc_all.m file calls helper functions which are 
contained in the helper_functions directory

Data_File_S1(2).txt:  Supplementary Files S1(2) from Richiardi et. al.

This code was updated in fall of 2019 in response the reply by Richiardi et. al.
to our 2017 commentary. The updates include new random network simulations with longer 
edge distance distributions and simulations that examine the effects of resting state
fMRI on significant SF frequencies among the simulated networks. The updated code
is mostly contained the last 3 cells of the proc_all.m file as well as in the 
helper_functions/random_clusters.m file. 

The Python code and Jupyter Notebook for running the spatial autocorrelation analyses are contained in the [spatial_autocorrelations](https://github.com/spiropan/ABA_functional_networks/blob/master/spatial_autocorrelations/spatial_autocorrelation.ipynb) subfolder. 

