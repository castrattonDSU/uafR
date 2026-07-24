# uafR offline student bundle verification.
# This script does not query live web services.

script_dir <- (function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    path <- sub("^--file=", "", file_arg[[length(file_arg)]])
    return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
  }
  frames <- sys.frames()
  for (frame in rev(frames)) {
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile, winslash = "/", mustWork = TRUE)))
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(path)) {
      return(dirname(normalizePath(path, winslash = "/", mustWork = TRUE)))
    }
  }
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
})()

helper_path <- file.path(script_dir, "bundle_helpers.R")
if (!file.exists(helper_path)) {
  helper_path <- file.path(script_dir, "..", "bundle_helpers.R")
}
source(helper_path)

bundle_dir <- uafr_bundle_root(script_dir)
log_file <- uafr_log_file(bundle_dir, "uafR_verify_log")
uafr_header("uafR Student Bundle Verification", log_file = log_file)

required_packages <- c("uafR", "ChemmineR", "fmcsR", "jsonlite", "webchem")
package_status <- uafr_package_status(required_packages)
uafr_print_status(package_status, log_file = log_file)

missing <- package_status$Package[!package_status$Installed]
if (length(missing) > 0) {
  uafr_note("Verification result: STOP. Missing package(s): ",
            paste(missing, collapse = ", "), log_file = log_file)
  stop("Verification failed. Missing package(s): ", paste(missing, collapse = ", "),
       call. = FALSE)
}

suppressPackageStartupMessages(library(uafR))

data_env <- new.env(parent = emptyenv())
datasets <- c("library_data", "standard_data", "standard_spread",
              "standard_exacto", "standard_categorated")
for (dataset in datasets) {
  utils::data(list = dataset, package = "uafR", envir = data_env)
  if (!exists(dataset, envir = data_env, inherits = FALSE)) {
    stop("Could not load package dataset: ", dataset, call. = FALSE)
  }
}

checks <- data.frame(
  Check = c(
    "R version",
    "Platform",
    "uafR version",
    "library_data rows",
    "library_data columns",
    "standard_data rows",
    "standard_spread entries",
    "standard_exacto rows",
    "standard_categorated tables"
  ),
  Value = c(
    R.version.string,
    R.version$platform,
    as.character(utils::packageVersion("uafR")),
    nrow(data_env$library_data),
    ncol(data_env$library_data),
    nrow(data_env$standard_data),
    length(data_env$standard_spread),
    nrow(data_env$standard_exacto),
    length(data_env$standard_categorated)
  ),
  stringsAsFactors = FALSE
)

uafr_note("", log_file = log_file)
print(checks, row.names = FALSE)
capture.output(print(checks, row.names = FALSE), file = log_file, append = TRUE)

research_schema_markers <- c("ValidationSummary", "TableQuality",
                             "SourceCoverage", "ChemicalTraitReport")
if (any(research_schema_markers %in% names(data_env$standard_categorated)) &&
    exists("validateCategorateResult", mode = "function")) {
  audit <- validateCategorateResult(data_env$standard_categorated)
  if (!is.null(audit$Summary)) {
    uafr_note("", log_file = log_file)
    print(audit$Summary)
    capture.output(print(audit$Summary), file = log_file, append = TRUE)
  }
} else {
  uafr_note("", log_file = log_file)
  uafr_note("Saved standard_categorated teaching data loaded successfully.",
            log_file = log_file)
  uafr_note("This object uses the legacy bundled teaching schema, so the research-schema validator is skipped.",
            log_file = log_file)
}

uafr_note("", log_file = log_file)
uafr_note("Verification result: PASS.", log_file = log_file)
uafr_note("The package loads, required dependencies are available, and bundled data are accessible.",
          log_file = log_file)
uafr_note("Log file: ", log_file, log_file = log_file)
