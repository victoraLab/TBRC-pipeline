# R packages installed from CRAN on top of the conda environment.
#
# environment.yml pins the conda builds (ggplate 0.1.3, ggplot2 3.5.1, ...). The production
# environment runs newer versions of ten packages, installed with install.packages() into the
# environment's own R library, which `conda env export` does not record. ggplate 0.1.3 fails in
# plate QC (approxfun, in its colour function), so run this once after `conda env create`.
#
#   conda activate snakemake && Rscript install_r_pins.R
#
# Activation matters: the packages are compiled with the environment's own compilers.
#
# Versions are those of the production environment on 2026-10-11.
pins <- c(ggplate = "0.2.0", ggplot2 = "4.0.2", scales = "1.4.0", dplyr = "1.2.0",
          readr = "2.2.0", vroom = "1.7.0", vctrs = "0.7.1", rlang = "1.1.7",
          gtable = "0.3.6", lifecycle = "1.0.5")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
for (p in names(pins)) {
  if (!identical(as.character(packageVersion(p)), pins[[p]]))
    remotes::install_version(p, version = pins[[p]], repos = "https://cloud.r-project.org", upgrade = "never")
}
chk <- vapply(names(pins), function(p) as.character(packageVersion(p)), "")
bad <- names(pins)[chk != pins]
if (length(bad)) stop("not at the pinned version: ", paste(bad, chk[bad], sep = " ", collapse = ", "), call. = FALSE)
cat("R pins in place:", paste(names(chk), chk, collapse = ", "), "\n")
