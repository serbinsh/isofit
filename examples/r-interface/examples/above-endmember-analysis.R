library(here)
library(raster)
library(tidyverse)
source(here("examples", "r-interface", "config.R"))
source(here("examples", "r-interface", "functions.R"))

reticulate::use_condaenv(CONDA_ENV)

envi_dir <- file.path("~", "projects", "sbg-uncertainty",
                      "sbg-prosail-workflow", "data",
                      "outputs", "envi-endmembers")

cwd <- here("examples", "r-interface", "examples")
atmosphere_df <- read.csv(file.path(cwd, "atmospheres.csv"))
outdir <- here("examples", "r-interface", "output")

endmembers <- c("Coastal_water", "Inland_water", "Mineral",
                "Snow", "Vegetation")

endmember_files <- file.path(envi_dir, endmembers)
wavelengths <- as.numeric(readLines(file.path(cwd, "above-wavelengths.txt")))
lrt <- read_libradtran_template(here(
  "examples", "r-interface", "lrt_template.inp"
))

# Loop over files
efile <- endmember_files[1]
dat <- t(brick(efile)[])
dat[dat < 0] <- 0

# Loop over spectra
outlist <- list()
for (i in seq_len(nrow(atmosphere_df))) {
  tag <- paste0("atmosphere", i)
  outlist[[tag]] <- list()
  for (r in seq_len(ncol(dat))) {
    reflectance <- dat[,r]
    outlist[[tag]][[r]] <- ht_workflow(
      reflectance,
      atmosphere_df[i, "H2O"],
      atmosphere_df[i, "AOD"],
      wavelengths,
      lrt,
      libradtran_basedir = LIBRADTRAN_DIR,
      outdir = outdir
    )
  }
}

reflmat <- function(x) do.call(cbind, lapply(x, "[[", "reflectance"))
results <- do.call(cbind, lapply(outlist, reflmat))

matplot(wavelengths, results, type = "l")
