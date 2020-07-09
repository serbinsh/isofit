# Isofit/Hypertrace R interface

## Installation

### Isofit

We recommend installing Isofit into a new conda environment.
Therefore, make sure you have anaconda or miniconda installed.
In this document, we'll call the conda environment `r-isofit`, but you can call it whatever you'd like (just replace `r-isofit` with your own environment name everywhere in these instructions).

First, download Isofit and checkout this branch.

``` sh
git clone https://github.com/ashiklom/isofit
cd isofit
git checkout r-hypertrace
```

Create a new conda environment and activate it.

``` sh
conda create -n r-isofit python
# Go through the interactive prompts...
conda activate r-isofit
```

Install Isofit and its dependencies.

``` sh
# From inside the isofit root directory (same directory as README.rst, LICENSE, etc.)
pip install .
```

### Libradtran

To generate your own atmospheric look-up tables, you'll need a working installation of LibRadTran.
Follow the instructions in the [Isofit `README`](https://github.com/ashiklom/isofit/tree/r-geom-2#quick-start-with-libradtran-20x) to install.
Note that you must install LibRadTran into the source code directory for it to work properly; i.e.

``` sh
# From the LibRadTran source code directory:
./configure --prefix=$(pwd)
make
```

### R packages

This code will require you to have the following R packages.
All of these are available on CRAN and can be installed via `install.packages` _except_ `rrtm` (optional, but used in the examples to run RTMs to simulate surface reflectance).

- `reticulate` -- For calling Python functions from R
- `R.matlab` -- Matlab IO (for writing the prior files)
- `digest` -- For creating hashes of the libradtran input file
- `fs` -- Better file path handling
- `here` -- Quickly locating files relative to a project root, regardless of working directory
- `rrtm` -- R interface to vegetation radiative transfer models (RTMs). Install from GitHub with `remotes::install_github("ashiklom/rrtm")` (will need to install `remotes` from CRAN if you don't already have it.)

## Configuration

In the `examples/r-interface` directory, copy the `config.example.R` file to `config.R` and modify accordingly.

## Running the workflow

For a complete, basic workflow, look at the `ht-workflow-fun.R` script, which generates a single surface reflectance spectrum and runs it through the "hypertrace" workflow:
i.e. simulate TOA radiance, and then atmospherically correct this radiance to estimate surface reflectance.
See the other example scripts for other workflows.
