# uafR private student bundle update script.
# Use this when a newer bundle is distributed.

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
log_file <- uafr_log_file(bundle_dir, "uafR_update_log")

uafr_header("uafR Student Bundle Update", log_file = log_file)
uafr_note("This script reinstalls uafR from this bundle and verifies the result.",
          log_file = log_file)
if (requireNamespace("uafR", quietly = TRUE)) {
  uafr_note("Current installed uafR version: ",
            as.character(utils::packageVersion("uafR")), log_file = log_file)
} else {
  uafr_note("uafR is not currently installed. This will run as a normal install.",
            log_file = log_file)
}

tryCatch(
  {
    uafr_install_from_bundle(bundle_dir, log_file = log_file, action = "update")
    uafr_note("", log_file = log_file)
    uafr_note("Update complete.", log_file = log_file)
    uafr_note("Restart RStudio if RStudio does not immediately find library(uafR).",
              log_file = log_file)
    uafr_note("Log file: ", log_file, log_file = log_file)
  },
  error = function(e) {
    uafr_note("", log_file = log_file)
    uafr_note("Update failed: ", conditionMessage(e), log_file = log_file)
    uafr_note("Log file: ", log_file, log_file = log_file)
    stop(e)
  }
)
