library(here)
source(here("examples/r-interface/zz-common.R"))
r <- ht_workflow(
    reflectance, 1.5, 0.2, wavelengths,
    libradtran_template,
    LIBRADTRAN_DIR,
    libradtran_environment = LIBRADTRAN_ENV,
    outdir = outdir
)

plot(wavelengths, r$reflectance, type = "l")
