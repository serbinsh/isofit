library(here)
library(raster)
source(here("examples", "r-interface", "config.R"))
source(here("examples", "r-interface", "functions.R"))

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

datlist <- lapply(endmember_files, function(x) t(brick(x)[]))
dat <- do.call(cbind, datlist)
dat[dat < 0] <- 0

outlist <- list()
for (i in seq_len(nrow(atmosphere_df))) {
  tag <- paste0("atmosphere", i)
  outlist[[tag]] <- ht_workflow(
      dat,
      atmosphere_df[i, "H2O"],
      atmosphere_df[i, "AOD"],
      wavelengths,
      lrt,
      libradtran_basedir = LIBRADTRAN_DIR,
      outdir = outdir
  )
}

results <- do.call(cbind, lapply(outlist, "[[", "reflectance"))
matplot(wavelengths, results, type = "l")
