library(here)
source(here("examples/r-interface/zz-common.R"))
# Default zenith angle
r1 <- ht_workflow(
    reflectance, 1.5, 0.2, wavelengths,
    libradtran_template,
    LIBRADTRAN_DIR,
    libradtran_environment = LIBRADTRAN_ENV,
    outdir = outdir
)

# Modify observer zenith angle
r2 <- ht_workflow(
    reflectance, 1.5, 0.2, wavelengths,
    libradtran_template,
    geom = list(observer_zenith = 20),
    LIBRADTRAN_DIR,
    libradtran_environment = LIBRADTRAN_ENV,
    outdir = outdir
)

matplot(wavelengths, cbind(r1$reflectance, r2$reflectance), type = "l")
