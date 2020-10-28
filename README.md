# OpenSystemZrSat.m

Matlab scripts to implement the stochastic fractional crystallization (FC) and recharge fraction crystallization (RFC) zircon saturation model from Tavazzani et al. (2020).

## Usage

The folder [LiteratureDatasets](LiteratureDatasets/) includes a compilation of CA-ID-TIMS zircon ages datasets for a number of intrusive and volcanic magmatic systems, stored as comma-separated files. A Matlab scrip to load and calculate kernel density estimates (KDEs) of all literature datasets is provided in [LiteratureDatasets/IDTIMS_DataLoad_KDE.m](LiteratureDatasets/IDTIMS_DataLoad_KDE.m). New datasets of interest can be added as comma-separated files in [LiteratureDatasets](LiteratureDatasets/), to load and and calculate KDEs of newly added zircon dates, the script [IDTIMS_DataLoad_KDE.m](LiteratureDatasets/IDTIMS_DataLoad_KDE.m) requires editing. 

For info and suggestions you can contact Lorenzo at [ltavazzani@smu.edu](mailto:ltavazzani@smu.edu).
