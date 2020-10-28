# OpenSystemZrSat.m

Matlab scripts to implement the stochastic fractional crystallization (FC) and recharge fraction crystallization (RFC) zircon saturation model from Tavazzani et al. (2020).

## Usage

The folder [SyntheticZirconDistribution](SyntheticZirconDistribution/) includes necessary scripts and data files necessary to obtain randomly sampled synthetic zircon age spectra from synthetic zircon crystallization distributions. 

In this folder, the files [dZr_Matlab_set2.xlsx](SyntheticZirconDistribution/dZr_Matlab_set2.xlsx) and [dZr_Matlab_set2_cutoff.xlsx](SyntheticZirconDistribution/dZr_Matlab_set2_cutoff.xlsx) contain zircon crystallization distributions for FC and RFC crystallization simulations, which are calculated offline from thermodynamic parameters output of the software [Magma Chamber Simulator](https://mcs.geol.ucsb.edu/code). 

To run the stochastic modeling and produce graphical outputs from several different crystallization distributions open the script [StochasticZirconCrystallization.m](SyntheticZirconDistribution/StochasticZirconCrystallization.m) and follow the directions. The non-uniform random generation of zircon age is performed by the function [StoZrnSamp.m](SyntheticZirconDistribution/StoZrnSamp.m), which does not require any editing by users.
 
The folder [LiteratureDatasets](LiteratureDatasets/) includes a compilation of CA-ID-TIMS zircon ages datasets for a number of intrusive and volcanic magmatic systems, stored as comma-separated files. A Matlab script to load and calculate kernel density estimates (KDEs) of all literature datasets is provided in [LiteratureDatasets/IDTIMS_DataLoad_KDE.m](LiteratureDatasets/IDTIMS_DataLoad_KDE.m). 

New datasets of interest can be added as comma-separated files in [LiteratureDatasets](LiteratureDatasets/). To load and and calculate KDEs of newly added zircon dates, the script [IDTIMS_DataLoad_KDE.m](LiteratureDatasets/IDTIMS_DataLoad_KDE.m) requires editing. 

For info and suggestions you can contact Lorenzo at [ltavazzani@smu.edu](mailto:ltavazzani@smu.edu).
