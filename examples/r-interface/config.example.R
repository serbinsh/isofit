# The conda environment where Isofit is installed
CONDA_ENV <- "r-isofit"
reticulate::use_condaenv(CONDA_ENV)

# The root directory of the LibRadTran source code
LIBRADTRAN_DIR <- "/path/to/libRadTran-2.0.3"

# Additional environment declarations necessary for LibRadTran to work May be
# necessary on some HPCs if you compiled LibRadTran against libraries in
# non-standard locations. Multiple statements should be separated by a newline
# character (`\n`).
LIBRADTRAN_ENV <- ""
# For example:
# LIBRADTRAN_ENV <- "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/.local/lib\nexport PATH=$PATH:~/.local/bin"
