# Improving inference in wastewater-based epidemiology by modelling the statistical features of digital PCR

Adrian Lison (1,2), Timothy R. Julian (3, 4, 5), and Tanja Stadler (1,2)

(1) ETH Zurich, Department of Biosystems Science and Engineering, Zurich, Switzerland\
(2) SIB Swiss Institute of Bioinformatics, Lausanne, Switzerland\
(3) Eawag, Swiss Federal Institute of Aquatic Science and Technology, 8600 Dubendorf, Switzerland\
(4) Swiss Tropical and Public Health Institute, 4123 Allschwil, Switzerland\
(5) University of Basel, 4055 Basel, Switzerland\
(*) Corresponding author: adrian.lison@bsse.ethz.ch

## Abstract
The growing field of wastewater-based infectious disease surveillance relies on
the quantification of pathogen concentrations in wastewater using polymerase
chain reaction (PCR) techniques. However, existing models for monitoring
pathogen spread using wastewater have often been adapted from methods for case
count data and neglect the statistical features of PCR techniques. In this
paper, we seek to overcome the widespread simplistic modelling of wastewater PCR
measurements as normally or log-normally distributed by proposing an appropriate
model for digital PCR (dPCR). Building on established statistical theory of
dPCR, we derive approximations for the coefficient of variation of measurements
and the probability of non-detection and propose a hurdle model-based likelihood
for estimating concentrations from dPCR measurements. Using simulations and
real-world data, we show that simple likelihoods based on normal or log-normal
distributions are misspecified, affecting the estimation of pathogen
concentrations and infection trends over time. In contrast, the proposed
dPCR-specific likelihood accurately models the distribution of dPCR measurements
and improves epidemiological estimates and forecasts even if details of the
laboratory protocol are unknown. The method has been implemented in the
open-source R package [EpiSewer](https://adrian-lison.github.io/EpiSewer/) to
improve wastewater-based monitoring of pathogens.

## Contents of this repository
*Code version: v1.0.0*

This repository contains the simulation code, data, and analysis scripts of the
study "Improving inference in wastewater-based epidemiology by modelling the
statistical features of digital PCR". The code can be used to reproduce all
figures and numerical results in the paper.

The repository is structured as follows:
- **code**: Contains utility functions and scripts.
- **data/ww_data**: Contains the real-world wastewater data used in this study.
  See the [data README](data/ww_data/README.md) for details.
- **data/results**: Contains pre-computed results of the simulation study. These
  can be used to reproduce the figures without re-running the simulations.
- **notebooks**: Contains all analysis scripts in R notebook format.
- **pipelines**: Contains scripts to fit wastewater models in `EpiSewer` using a
  modeling pipeline via the `targets` package.
- **renv**: Contains the configuration of `renv` for the project.

### Setup
The analysis scripts are written as R notebooks and are ideally run in an Rstudio
project. When opening the project, run `renv::restore()` to install all required
R packages (requires the `renv` package).

Note that [EpiSewer v0.0.3](https://doi.org/10.5281/zenodo.13899759) was used
for inference from concentration measurements in this study. To ensure
reproducibility, `renv` will automatically install this version of the package.
Newer versions of the package may be found on the [EpiSewer GitHub
repository](https://github.com/adrian-lison/EpiSewer).

### Simulation study
To reproduce the model validation via simulation, run the notebook
[CV_p_nondetect_sim.Rmd](notebooks/dPCR%20paper/CV_p_nondetect_sim.Rmd).

### Real-world validation
To reproduce the comparison with empirical data, run the notebook 
[CV_p_nondetect_real_world_resp6.Rmd](notebooks/dPCR%20paper/CV_p_nondetect_real_world_resp6.Rmd).

### Inference from concentration measurements
*Note: To rerun the model fitting, the stan inference engine must be installed.
See [the cmdstanr vignette](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)
for help.*

To reproduce the inference from measurements of a single concentration, run the
notebook [single_concentration.Rmd](notebooks/dPCR%20paper/single_concentration.Rmd).

To reproduce the inference of concentrations over time, run the notebook
[noise_model_comparison.Rmd](notebooks/dPCR%20paper/noise_model_comparison.Rmd).

### Additional results

- [averaging_measurements.Rmd](notebooks/dPCR%20paper/averaging_measurements.Rmd) 
  shows a comparison of concentration estimates from multiple replicates obtained 
  by simple averaging vs. by computing the ML estimate from the combined partition counts.
- [likelihoods_continuous.Rmd](notebooks/dPCR%20paper/likelihoods_continuous.Rmd) 
  compares the suitability of different continuous distributions for the
  likelihood function of the dPCR model.
- [nonzero_conditioning.Rmd](notebooks/dPCR%20paper/nonzero_conditioning.Rmd) 
  shows the effect of conditioning on non-zero measurements.
- [partition_loss.Rmd](notebooks/dPCR%20paper/partition_loss.Rmd) shows the 
  modeling of partition loss.

### Useful scripts
The script [utils_dPCR_statistics.R](code/utils_dPCR_statistics.R) contains
various functions to compute the coefficient of variation and the probability of
non-detection for dPCR measurements.

The script [_noise_model_comparison.R](pipelines/_noise_model_comparison.R) and
the utility scripts in [code/pipeline](code/pipeline) show how to fit wastewater
models in [EpiSewer](https://adrian-lison.github.io/EpiSewer/) using a modeling
pipeline via the [targets](https://books.ropensci.org/targets/) package.