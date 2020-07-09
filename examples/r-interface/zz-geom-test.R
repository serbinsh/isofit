library(here)
source(here("examples/r-interface/zz-common.R"))
source(here("examples/r-interface/functions.R"))
r <- ht_workflow(
    reflectance, 1.5, 0.2, wavelengths,
    libradtran_template,
    libradtran_basedir,
    lut = lut,
    outdir = outdir
)

plot(wavelengths, r$reflectance, type = "l")
