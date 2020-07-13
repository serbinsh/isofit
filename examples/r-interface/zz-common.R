library(reticulate)
library(here)
library(fs)
library(digest)

source(here("examples", "r-interface", "config.R"))
source(here("examples", "r-interface", "functions.R"))

use_condaenv(CONDA_ENV)
template_file <- here("examples", "r-interface", "lrt_template.inp")

# Arguments
wavelengths <- seq(400, 2500, 5)
true_refl <- drop(rrtm::pro4sail_4(1.4, 40, 0.01, 0.01, 3, 0.5)$bdr)
reflectance <- approx(400:2500, true_refl, wavelengths,
                      yleft = 0, yright = 0)$y
libradtran_template <- read_libradtran_template(template_file)

# Outputs stored in this directory.
outdir <- dir_create(here("examples", "r-interface", "output"))
