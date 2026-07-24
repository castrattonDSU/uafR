message("WiCCED water-quality training setup check")
message("R version: ", R.version.string)
message("Working directory: ", getwd())

required <- c("stats", "utils", "graphics")
optional <- c("uafR", "testthat", "ggplot2", "dplyr", "ranger", "randomForest")

check_pkg <- function(pkg) {
  available <- requireNamespace(pkg, quietly = TRUE)
  version <- if (available) as.character(utils::packageVersion(pkg)) else NA_character_
  data.frame(package = pkg, available = available, version = version,
             stringsAsFactors = FALSE)
}

status <- do.call(rbind, lapply(c(required, optional), check_pkg))
print(status, row.names = FALSE)

if (!requireNamespace("uafR", quietly = TRUE)) {
  message("uafR is not currently available to this R session.")
  message("Use the local student bundle or repository installation instructions before running uafR labs.")
} else {
  message("uafR detected.")
}

utils::write.csv(status, "setup_package_status.csv", row.names = FALSE)
message("Wrote setup_package_status.csv")
