# ParesPujolras_Kelly_Murphy_2024

This repository contains analysis scripts for the following paper:

[Parés-Pujolràs, E., Kelly, S. P., & Murphy, P. R. (2025). Dissociable encoding of evolving beliefs and momentary belief updates in distinct neural decision signals. Nature Communications, 16(1), 3922.]

The scripts included in the repository and the data on OSF (https://doi.org/10.17605/OSF.IO/NHQEV) enable the replication of the data analysis and figures reported in the paper.
All are shared under a CC-BY 4.0 license.

Code shared here was developed using Matlab R2021a, the [EEGlab toolbox](https://sccn.ucsd.edu/eeglab/index.php) (v.2021.0), [FieldTrip toolbox](https://www.fieldtriptoolbox.org/), [CSD toolbox](https://psychophysiology.cpmc.columbia.edu/software/csdtoolbox/), and the [TBT toolbox](https://github.com/mattansb/TBT).

## Behavioural analysis:
Contains scripts for fitting several models to choice data & extracting transformed stimulus values based on the fitted models. 
The basic model-fitting functions are based on analysis code from a previous study ([Murphy et al. 2021, Nat. Neuro](https://www.nature.com/articles/s41593-021-00839-z)), with publicly available code [here](https://github.com/DonnerLab/2021_Murphy_Adaptive-Circuit-Dynamics-Across-Human-Cortex). 
Includes plotting code for all behavioural data & model fitting results in Figure 1.

## EEG analysis: 
Includes raw EEG data loading, preprocessing, regression & residuals analysis, as well as statistics and plotting code for Figures 2-4 in the paper.  


