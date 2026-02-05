# *Nc*AA9D Dose Series Analysis Scripts

The scripts and various files in this repository were used to process data and
perform statistical analyses for the article
["Dose-dependent structural and electron density features in the lytic polysaccharide monooxygenase NcAA9D"](LINK HERE).
This repository is organized with the intent of being transparent in methodology
and processing flow. Five directories—numbered by their order in the pipeline—hold
the scripts and majority of input data (file size allowing) used in our analyses.
If desired, following the instructions below should allow you to fully reproduce
the data used in our article.

Due to the large dataset size, our analysis relied on heavy use of Bash
and R scripts to organize files, improve efficiency, minimize human error,
and maximize reproducibility/transparency. Linux and Mac users may simply run
`bash path/to/script/script_name.sh` from the terminal (ensure your current working directory
is NcAA9D_DoseSeries-main). Windows users may run these scripts using a Bash
interpreter such as Git Bash. R may be downloaded from the (R Project website)[https://www.r-project.org/]
or from the (RStudio website)[https://posit.co/download/rstudio-desktop/].

If you would like to replicate any results obtained from R scripts, please open
r_analysis.Rproj and run the command `renv::restore()` (if you haven't used R before
you will need to download renv with `install.packages("renv")`) . The desired
script may then be run by typing the command `source("path/to/script/script_name.R")`
into the R terminal.

### 1. Diffraction Decay Parameter Estimation

Our unique data collection strategy required us to calculate the average
**diffraction decay-weighted dose** for each dataset. To calculate this in
RADDOSE-3D, we three parameters must be estimated:

1. γ: describes the dose-dependent behavior of the Gaussian scale factor (MGy<sup>-1</sup>)

2. B<sub>0</sub>: the Wilson B-factor at zero dose (Å<sup>2</sup>)

3. β: the rate of Wilson B-factor increase per unit dose (Å<sup>2</sup>/MGy)

The first 25 "sliding window" pseudohelices, corresponding to start angles
0-24° (0-5°, 1-6°, ... 24-29°) were processed in [DIALS](https://dials.github.io/index.html)
(v3.26). The integrated (unmerged and unscaled) reflection files were exported,
and their absolute scale factors and Wilson B-factors were estimated with
[WILSON](https://www.ccp4.ac.uk/html/wilson.html) (v9.0.011). Diffraction
decay parameters were then estimated by least-squares regression in R, according
to the method described in [Dickerson *et al.* (2024)](https://doi.org/10.1002/pro.5005).
Further details of this analysis are available in our article.

To run the analysis from scratch, you will need the diffraction images. These are available
at **LINK HERE**, and should be extracted into 0-diffraction_images. You may then
run the data processing script at 1-decay_parameter_estimation/...........sh
This will process the necessary datasets and output a CSV file with the datasets'
Wilson B-factors and scale factors.

If you don't have the computer storage or time to spare, no worries! All the
input files you need to estimate γ, B<sub>0</sub>, and β are already here
(1-decay_parameter_estimation/input/wilson_b_scales.csv). The R script
to do this is located at 1-decay_parameter_estimation/estimate_gamma_b0_beta.R.
The script will take the input Wilson B-factors and scale factors,
estimate γ, B<sub>0</sub>, and β, export the values as a CSV file, and export
the relevant plots.

### 2. Dose Calculations

After γ, B<sub>0</sub>, and β were estimated, these were used to calculate the
average diffraction decay-weighted doses (DDWDs) for our wedge datasets in
[RADDOSE-3D](https://github.com/GarmanGroup/RADDOSE-3D) (v5.01058). The average 
fluence-weighted doses (FWDs), which don't require these parameters, were also
calculated. Finally, an R script was used to calculate the FWDs and DDWDs for
the pseudohelix datasets, as this can't be done natively with RADDOSE.

To calculate DDWDs and FWDs for all the wedge datasets, run 2-dose_calculations/..........sh
This will generate the relevant RADDOSE-3D input files, calculate the wedge FWDs
and DDWDs with RADDOSE-3D, and export the relevant output files.

All the input files necessary to calculate pseudohelix doses are already here!
Simply run 2-dose_calculations/fwd_ddwd_calculation.R, which will calculate average
FWDs/DDWDs and export a CSV file and plot of all the results.

### 3. Sampling

Weighted random sampling was used to select a set of datasets to use in 
analysis. Higher weights were given to pseudohelices with higher rates of 
average DDWD accumulation, with the goal of obtaining a sample set
representative of the whole dose range achieved in our study.

The input files you need to reproduce this step were originally produced by the
scripts in section 2, but are all here for your convenience. To reproduce the
sampling procedure, simply run 3-sampling/sampling.R. This script performs the
weighted random sampling procedure and exports a CSV file and plot of the results.

Note that sampling.R hard-codes the random seed to produce the same sample set
used in our article. If you'd like to see the random sampling in action, remove
the line `set.seed(12)` from the script and run it again.

### 4. Data Processing



### 5. Structure Refinements

After processing the desired datasets, the relevant reflections files and a
"base" *Nc*AA9D model were used as an input in [phenix.refine](https://phenix-online.org/documentation/reference/refinement.html)
(vXXXXXXXX) to determine the structure corresponding to each dataset. This is a computationally
hefty process, so we don't recommend you run the relevant scripts. We have, though,
included the configuration files, Bash scripts, and file structure used in this
step of the process to improve clarity.

### 6. Structure Analysis

The analysis pipeline built for this analysis imports and parses the .pdb files, 
storing structure parameters in a large table. Multiple linear regression models
are then calculated for parameters of interest versus dose and chain ID (see
the publication for more details). The script also generates plots which were
used in the paper.
