library(here)
library(raster)
source(here("examples", "r-interface", "config.R"))
source(here("examples", "r-interface", "functions.R"))

cwd <- here("examples", "r-interface", "examples")
atmosphere_df <- read.csv(file.path(cwd, "atmospheres.csv"))
outdir <- here("examples", "r-interface", "output")
lrt <- read_libradtran_template(here(
  "examples", "r-interface", "lrt_template.inp"
))

endmembers <- read.csv(file.path(cwd, "endmembers.csv"))
wavelengths <- endmembers$wavelength
endmember_mat <- as.matrix(endmembers[,-1])

outlist <- list()
for (i in seq_len(nrow(atmosphere_df))) {
  tag <- paste0("atmosphere", i)
  outlist[[tag]] <- ht_workflow(
    endmember_mat,
    atmosphere_df[i, "H2O"],
    atmosphere_df[i, "AOD"],
    wavelengths,
    lrt,
    outdir
  )
}

results <- do.call(cbind, lapply(outlist, "[[", "reflectance"))
matplot(wavelengths, results, type = "l")
