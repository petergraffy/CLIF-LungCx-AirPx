# 00_install_packages.R
# Install all required packages for this project.
#
# Run this script ONCE before running any other scripts.
# Only installs what's missing — safe to re-run.

required_packages <- c(
  # Core tidyverse + data manipulation
  "tidyverse", "stringr", "lubridate", "readr", "glue",
  "dplyr", "data.table", "scales",
  # File I/O + cleaning
  "arrow", "fst", "jsonlite", "janitor",
  # Clustering & trajectory analysis (02_federated_clusters.R)
  "TraMineR", "cluster", "nnet", "MASS", "survival", "comorbidity",
  # Spatial / mapping
  "sf", "tigris"
)

missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  install.packages(missing)
} else {
  message("All required packages are already installed.")
}

# Verify
check <- vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)
if (all(check)) {
  message("\n=== Setup complete. You can now run 01_lungcx_cohort.R ===")
} else {
  warning(
    "Failed to install: ",
    paste(required_packages[!check], collapse = ", "),
    "\nPlease install them manually."
  )
}
