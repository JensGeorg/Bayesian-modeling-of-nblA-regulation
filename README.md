# Bayesian-modeling-of-nblA-regulation

This repository contains a Stan model for fitting a system of ordinary differential equations (ODEs) that describe the dynamics of the *nblA* (RNA), as_nblA (asRNA), and NsrR1 (sRNA) regulatory system. The model is fitted to various experimental datasets, including microarray time-series data, Northern blot data, and RNA decay rates. The repository also includes R scripts to process the fitted model, visualize the results, and run new simulations for different biological scenarios.

## Model Description

The Stan model `rna_dynamics.stan` implements a four-state ODE system:
1.  `RNA_free`: Unbound nblA concentration.
2.  `RNA_complex_sRNA`: nblA/NsrR1 complex concentration.
3.  `asRNA_free`: Unbound as_nblA concentration.
4.  `sRNA_free`: Unbound NsrR1 concentration.

## Repository Contents

* `stan/rna_dynamics.stan`: The core Stan code for the Bayesian model.
* `R/plot_fits_to_data.r`: An R script that loads the fitted Stan object and generates a plot comparing the model's posterior predictions with the observed data.
* `R/extract_fitted_pars_tocsv.r`: An R script to extract the posterior distributions of all model parameters and save them to a CSV file.
* `R/run_simulations.r`: An R script to run forward simulations of the ODE model using the mean posterior parameter estimates. It simulates different biological scenarios, such as wild type and knockout conditions.
* `data/`: Directory to store any raw data files (not included in this template, but recommended).
* `run_analysis.R`: A master script to load the data, fit the model, and run the analysis scripts.
* `.gitignore`: A file to specify files to be ignored by Git (e.g., large output files, `.RData` files).

## Getting Started

### Prerequisites

* R environment
* Stan (`rstan` package)
* Required R libraries: `ggplot2`, `dplyr`, `tidyr`, `deSolve`, `boot`, `ggsci`

### Installation

1.  Clone this repository to your local machine.
    ```sh
    git clone [https://github.com/your-username/your-repo-name.git](https://github.com/your-username/your-repo-name.git)
    cd your-repo-name
    ```
2.  Install the necessary R packages.
    ```R
    install.packages(c("rstan", "ggplot2", "dplyr", "deSolve", "tidyr", "boot", "ggsci"))
    ```

### Usage

1.  Place your data files in the `data/` directory (if applicable).
2.  Modify the data loading and list creation in the `run_analysis.R` script if your data format is different.
3.  Run the main analysis script from the R console. This script will perform the model fitting and subsequent analysis.
    ```R
    source("run_analysis.R")
    ```

The `run_analysis.R` script is configured to:
-   Define all data for the Stan model.
-   Compile and fit the `rna_dynamics.stan` model.
-   Save the fitted parameters to `fitted_parameters_and_priors.csv`.
-   Generate `timeseries.pdf` showing the model fit to data.
-   Generate `test35.pdf` showing the simulated scenarios.

## Scripts

### `stan/rna_dynamics.stan`

(The cleaned Stan code you provided goes here.)

### `R/plot_fits_to_data.r`

This script is for visualizing the model fit. It expects a fitted Stan object named `fit` to be available in the R environment.

```R
# R code for plot_fits_to_data.r goes here
