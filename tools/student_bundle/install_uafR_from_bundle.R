# uafR private student bundle installer.
# Unzip the bundle first, then open this file in RStudio and click Source.

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
log_file <- uafr_log_file(bundle_dir, "uafR_install_log")

uafr_header("uafR Student Bundle Install", log_file = log_file)
uafr_note("This installer uses CRAN and Bioconductor for dependencies.",
          log_file = log_file)
uafr_note("The uafR package itself is installed from this local bundle.",
          log_file = log_file)

tryCatch(
  {
    uafr_install_from_bundle(bundle_dir, log_file = log_file, action = "install")
    uafr_note("", log_file = log_file)
    uafr_note("Install complete.", log_file = log_file)
    uafr_note("Restart RStudio if RStudio does not immediately find library(uafR).",
              log_file = log_file)
    uafr_note("Log file: ", log_file, log_file = log_file)
  },
  error = function(e) {
    uafr_note("", log_file = log_file)
    uafr_note("Install failed: ", conditionMessage(e), log_file = log_file)
    uafr_note("Log file: ", log_file, log_file = log_file)
    stop(e)
  }
)
