# Improving inference in environmental surveillance by modelling the statistical features of digital PCR

Adrian Lison (1,2), Timothy R. Julian (3, 4, 5), and Tanja Stadler (1,2)

(1) ETH Zurich, Department of Biosystems Science and Engineering, Zurich, Switzerland\
(2) SIB Swiss Institute of Bioinformatics, Lausanne, Switzerland\
(3) Eawag, Swiss Federal Institute of Aquatic Science and Technology, 8600 Dubendorf, Switzerland\
(4) Swiss Tropical and Public Health Institute, 4123 Allschwil, Switzerland\
(5) University of Basel, 4055 Basel, Switzerland\
(*) Corresponding author: adrian.lison@bsse.ethz.ch

## Abstract
Digital polymerase chain reaction (dPCR) is a powerful technique for quantifying
gene targets in environmental samples, with various applications such as
biodiversity monitoring and wastewater-based epidemiology. However, statistical
analyses of environmental dPCR data often assume, explicitly or implicitly, that
concentration measurements have a normal or log-normal error structure, which
does not reflect the underlying partitioning statistics of dPCR. Using
simulations and real-world environmental data, we show that (log-)normality
assumptions are violated for dPCR measurements, leading to inaccurate estimates
of gene concentrations and underlying biological processes. To enable reliable
analyses of environmental dPCR data, we present a dPCR-specific likelihood model
that accounts for concentration-dependent measurement noise and non-detects as
characteristic of dPCR assays. We demonstrate that this approach overcomes
biases in inference from environmental data, such as estimating free-eDNA decay
in seawater or pathogen transmission from wastewater monitoring. Our method is
implemented in the R packages
[dPCRfit](https://github.com/adrian-lison/dPCRfit/) for regression analyses and
[EpiSewer](https://adrian-lison.github.io/EpiSewer/) for wastewater
surveillance.

## Contents of this repository
*Code version: v2.0.0*

This repository contains the simulation code, data, and analysis scripts of the
study "Improving inference in environmental surveillance by modelling the
statistical features of digital PCR". The code can be used to reproduce all
figures and numerical results in the paper.

The repository is structured as follows:
- **code**: Contains utility functions and scripts.
- **data/ww_data**: Contains real-world wastewater data used in this study.
  See the [wastewater data README](data/ww_data/README.md) for details.
- **data/eDNA**: Contains eDNA measurements by Scriver et al. that were reanalyzed in this study.
  See the [eDNA data README](data/eDNA/README.md) for details.
- **data/results**: Contains pre-computed results and fitted models. These
  can be used to reproduce the figures without re-running simulations / model fits.
- **notebooks**: Contains all analysis scripts in R notebook format.
- **pipelines**: Contains scripts to fit wastewater models in `EpiSewer` using a
  modeling pipeline via the `targets` package.
- **renv**: Contains the configuration of `renv` for the project.

### Setup
The analysis scripts are written as R notebooks and are ideally run in an Rstudio
project. When opening the project, run `renv::restore()` to install all required
R packages (requires the `renv` package).

Note that the R packages [dPCRfit v0.0.2](https://doi.org/10.5281/zenodo.15089580) 
and [EpiSewer v0.0.3](https://doi.org/10.5281/zenodo.13899759) were used
for inference from concentration measurements in this study. To ensure
reproducibility, `renv` will automatically install this version of the package.
Newer versions of the package may be found on the 
[dPCRfit](https://github.com/adrian-lison/dPCRfit) and 
[EpiSewer](https://github.com/adrian-lison/EpiSewer) Github repositories.

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

To reproduce the reanalysis of eDNA measurements by Scriver et al., run the
notebook [aquarium_Scriver_et_al.Rmd](notebooks/dPCR%20paper/aquarium_Scriver_et_al.Rmd).

To reproduce the analysis of wastewater measurements with different noise models, run the notebook
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