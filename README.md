# *Nc*AA9D Dose Series Analysis Scripts

The scripts and various files were used to process data and perform statistical
analysis for the article
["Dose-dependent structural and electron density features in the lytic polysaccharide monooxygenase NcAA9D"](LINK HERE).
This repository is organized with the intent of being transparent in methodology
and processing flow. Five directories hold the scripts and majority of input
data used in the analyses. If desired, following the instructions below should
allow you to fully reproduce the data used in the above article.

### 1. Diffraction Decay Parameter Estimation

Our unique data collection strategy required us to calculate the average
**diffraction decay-weighted doses** for each dataset. To calculate this in
RADDOSE-3D, we needed to estimate three parameters:

1. γ: describes the dose-dependent behavior of the Gaussian scale factor (MGy<sup>-1</sup>)

2. B<sub>0</sub>: the Wilson B-factor at zero dose (Å<sup>2</sup>)

3. β: the rate of Wilson B-factor increase per unit dose (Å<sup>2</sup>/MGy)

The first 25 "sliding window" pseudohelices, corresponding to start angles
0-24° (0-5°, 1-6°, ... 24-29°) were processed in DIALS (v3.26). The integrated
(unmerged, unscaled) reflection files were exported, and their absolute scale
factors and Wilson B-factors were calculated with WILSON (v9.0.011). The above
parameters were then estimated by least-squares regression, according to the
method described in [Dickerson *et al.* (2024)](https://doi.org/10.1002/pro.5005).
Further details of the model and parameters are available in our article.

INSERT PLOT HERE

To run the analysis, you will need the diffraction images. These are available
at LINK HERE, and should be extracted into `0-diffraction_images`. You may then
run the data processing script at `./1-diffration_decay_parameter_estimation/..........`

### 2. Dose Calculations



### 3. Sampling



### 4. Structure Refinements



### 5. Structure Analysis

The analysis pipeline built for this analysis imports and parses the .pdb files, 
storing structure parameters in a large table. Multiple linear regression models
are then calculated for parameters of interest versus dose and chain ID (see
the publication for more details). The script also generates plots which were
used in the paper.



The scripts are packed together as an R project. To run the analysis scripts, 
simply run renv::restore() upon opening Analysis.Rproj.
