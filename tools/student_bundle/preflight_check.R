# uafR private student bundle preflight check.
# Run this before installation if a machine has not used uafR before.

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
log_file <- uafr_log_file(bundle_dir, "uafR_preflight_log")
uafr_header("uafR Student Bundle Preflight Check", log_file = log_file)

rows <- list()
add_check <- function(check, status, details) {
  rows[[length(rows) + 1]] <<- uafr_status_row(check, status, details)
}

archive <- tryCatch(uafr_find_package_archive(bundle_dir), error = function(e) NA_character_)
add_check(
  "Bundle archive",
  if (!is.na(archive)) "PASS" else "FAIL",
  if (!is.na(archive)) basename(archive) else "No packages/uafR_<version>.tar.gz file was found."
)

required_files <- c(
  "install_uafR_from_bundle.R",
  "update_uafR_from_bundle.R",
  "verify_uafR_install.R",
  "run_student_acceptance_test.R",
  "README_STUDENT_INSTALL.md"
)
missing_files <- required_files[!file.exists(file.path(bundle_dir, required_files))]
add_check(
  "Bundle scripts",
  if (length(missing_files) == 0) "PASS" else "FAIL",
  if (length(missing_files) == 0) "Required bundle scripts are present." else
    paste("Missing:", paste(missing_files, collapse = ", "))
)

r_ok <- getRversion() >= "4.2.0"
add_check(
  "R version",
  if (r_ok) "PASS" else "FAIL",
  paste0(R.version.string, if (r_ok) "" else " -- update R before installing uafR.")
)

rstudio_ok <- requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()
add_check(
  "RStudio",
  if (rstudio_ok) "PASS" else "WARN",
  if (rstudio_ok) "RStudio is available." else
    "RStudio was not detected. The scripts can still run, but student instructions assume RStudio."
)

library_ok <- uafr_library_writable()
add_check(
  "Writable R library",
  if (library_ok) "PASS" else "FAIL",
  paste0(.libPaths()[1], if (library_ok) "" else " is not writable.")
)

cran_ok <- uafr_download_ok("https://cloud.r-project.org/src/contrib/PACKAGES.gz")
add_check(
  "CRAN access",
  if (cran_ok) "PASS" else "WARN",
  if (cran_ok) "Reached https://cloud.r-project.org." else
    "Could not reach CRAN. Installation can still work only if CRAN dependencies are already installed."
)

bioc_ok <- uafr_download_ok("https://bioconductor.org/config.yaml")
add_check(
  "Bioconductor access",
  if (bioc_ok) "PASS" else "WARN",
  if (bioc_ok) "Reached https://bioconductor.org." else
    "Could not reach Bioconductor. Installation can still work only if Bioconductor dependencies are already installed."
)

dependency_status <- uafr_package_status(c("ChemmineR", "fmcsR", "jsonlite", "webchem"))
missing_dependencies <- dependency_status$Package[!dependency_status$Installed]
add_check(
  "Current dependencies",
  if (length(missing_dependencies) == 0) "PASS" else "WARN",
  if (length(missing_dependencies) == 0) "All required dependencies are already installed." else
    paste("Missing now:", paste(missing_dependencies, collapse = ", "))
)

if (.Platform$OS.type == "windows") {
  make_path <- Sys.which("make")
  add_check(
    "Windows build tools",
    if (nzchar(make_path)) "PASS" else "WARN",
    if (nzchar(make_path)) paste("make found at", make_path) else
      "Rtools was not detected. Install it only if R says a package must be compiled from source."
  )
} else if (identical(Sys.info()[["sysname"]], "Darwin")) {
  make_path <- Sys.which("make")
  add_check(
    "macOS build tools",
    if (nzchar(make_path)) "PASS" else "WARN",
    if (nzchar(make_path)) paste("make found at", make_path) else
      "Command line build tools were not detected. Install them only if R says a package must be compiled from source."
  )
}

status_table <- do.call(rbind, rows)
uafr_note("", log_file = log_file)
uafr_print_status(status_table, log_file = log_file)

uafr_note("", log_file = log_file)
uafr_note("Dependency details:", log_file = log_file)
uafr_print_status(dependency_status, log_file = log_file)

failures <- status_table$Check[status_table$Status == "FAIL"]
warnings <- status_table$Check[status_table$Status == "WARN"]

uafr_note("", log_file = log_file)
if (length(failures) > 0) {
  uafr_note("Preflight result: STOP. Fix the failed item(s): ",
            paste(failures, collapse = ", "), log_file = log_file)
  stop("Preflight found blocking setup problems. See ", log_file, call. = FALSE)
}

if (length(warnings) > 0) {
  uafr_note("Preflight result: READY WITH WARNINGS. The installer can be tried, but review: ",
            paste(warnings, collapse = ", "), log_file = log_file)
} else {
  uafr_note("Preflight result: READY. This machine is ready to install uafR.",
            log_file = log_file)
}

uafr_note("Log file: ", log_file, log_file = log_file)
