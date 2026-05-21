# DSU dsDNA Core training install check.
# Run from the repository root or from a student project that has uafR installed.

check_training_install <- function() {
  required <- c("uafR", "ChemmineR", "fmcsR", "jsonlite", "webchem")

  package_status <- data.frame(
    Package = required,
    Installed = vapply(required, requireNamespace, logical(1), quietly = TRUE),
    Version = vapply(required, function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
      as.character(utils::packageVersion(pkg))
    }, character(1)),
    stringsAsFactors = FALSE
  )

  print(package_status, row.names = FALSE)

  if (!requireNamespace("uafR", quietly = TRUE)) {
    stop("uafR is not installed. Install uafR before continuing.", call. = FALSE)
  }

  library(uafR)

  data("library_data", package = "uafR")
  data("standard_data", package = "uafR")

  checks <- list(
    RVersion = R.version.string,
    Platform = R.version$platform,
    UafRVersion = as.character(utils::packageVersion("uafR")),
    LibraryDataRows = nrow(library_data),
    LibraryDataColumns = ncol(library_data),
    StandardDataRows = nrow(standard_data),
    StandardDataColumns = ncol(standard_data)
  )

  print(checks)
  invisible(list(packages = package_status, checks = checks))
}
