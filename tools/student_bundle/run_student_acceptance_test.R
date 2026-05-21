# uafR student acceptance test.
# By default this script avoids live web-service queries.
# To include the live database smoke test, run with --live or set UAFR_STUDENT_LIVE=1.

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
log_file <- uafr_log_file(bundle_dir, "uafR_acceptance_test_log")
args <- commandArgs(trailingOnly = TRUE)
live_requested <- any(args %in% c("--live", "live")) ||
  identical(tolower(Sys.getenv("UAFR_STUDENT_LIVE")), "1") ||
  identical(tolower(Sys.getenv("UAFR_STUDENT_LIVE")), "true")

uafr_header("uafR Student Acceptance Test", log_file = log_file)
uafr_note("Live database smoke test: ", if (live_requested) "ON" else "OFF",
          log_file = log_file)

results <- list()
run_check <- function(name, expr) {
  uafr_note("", log_file = log_file)
  uafr_note("Checking: ", name, log_file = log_file)
  detail <- tryCatch(
    {
      value <- force(expr)
      if (is.null(value)) "completed" else paste(value, collapse = "; ")
    },
    error = function(e) {
      structure(conditionMessage(e), class = "uafr_check_error")
    }
  )

  if (inherits(detail, "uafr_check_error")) {
    uafr_note("Result: FAIL - ", as.character(detail), log_file = log_file)
    results[[length(results) + 1]] <<-
      uafr_status_row(name, "FAIL", as.character(detail))
  } else {
    uafr_note("Result: PASS - ", detail, log_file = log_file)
    results[[length(results) + 1]] <<-
      uafr_status_row(name, "PASS", detail)
  }
  invisible(TRUE)
}

run_check("Package verification", {
  uafr_run_verification(bundle_dir, log_file = log_file)
  "package, dependencies, and bundled data verified"
})

run_check("Offline training workflow", {
  workflow_script <- file.path(bundle_dir, "training", "scripts", "04_core_workflow.R")
  if (!file.exists(workflow_script)) {
    stop("Training workflow script is missing: ", workflow_script, call. = FALSE)
  }
  source(workflow_script, local = TRUE)
  result <- run_core_workflow(live_lookup = FALSE)
  if (!is.list(result) || !"exact" %in% names(result)) {
    stop("Offline workflow did not return the expected result list.", call. = FALSE)
  }
  if (nrow(result$exact) == 0) {
    stop("Offline workflow returned no exact-match rows.", call. = FALSE)
  }
  paste0("offline exact-match rows: ", nrow(result$exact))
})

run_check("Core package data objects", {
  suppressPackageStartupMessages(library(uafR))
  data_env <- new.env(parent = emptyenv())
  utils::data("library_data", package = "uafR", envir = data_env)
  utils::data("standard_data", package = "uafR", envir = data_env)
  utils::data("standard_spread", package = "uafR", envir = data_env)
  utils::data("standard_exacto", package = "uafR", envir = data_env)
  stopifnot(is.data.frame(data_env$library_data))
  stopifnot(is.data.frame(data_env$standard_data))
  stopifnot(is.list(data_env$standard_spread))
  stopifnot(is.data.frame(data_env$standard_exacto))
  paste0("library rows: ", nrow(data_env$library_data),
         "; standard data rows: ", nrow(data_env$standard_data),
         "; exact rows: ", nrow(data_env$standard_exacto))
})

if (live_requested) {
  run_check("Optional live categorate smoke test", {
    suppressPackageStartupMessages(library(uafR))
    data("library_data", package = "uafR")
    live_warnings <- character()
    quick_result <- withCallingHandlers(
      categorate(
        compounds = c("aspirin", "caffeine"),
        chemical_library = library_data,
        input_format = "wide",
        detail = "research",
        cache = TRUE,
        throttle = 0.2,
        assay_detail_limit = 0
      ),
      warning = function(w) {
        live_warnings <<- c(live_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    audit <- validateCategorateResult(quick_result)
    if (!is.null(audit$Summary)) {
      capture.output(print(audit$Summary), file = log_file, append = TRUE)
    }
    required_tables <- c("ChemicalTraitReport", "ChemicalMeasurementSummary")
    missing_tables <- setdiff(required_tables, names(quick_result))
    if (length(missing_tables) > 0) {
      stop("Live categorate result is missing table(s): ",
           paste(missing_tables, collapse = ", "), call. = FALSE)
    }
    if (length(live_warnings) > 0) {
      uafr_note("Non-blocking live lookup warning(s): ", length(live_warnings),
                log_file = log_file)
      writeLines(unique(live_warnings), con = log_file, sep = "\n")
    }
    paste0("live result tables: ", length(quick_result),
           "; non-blocking warnings: ", length(unique(live_warnings)))
  })
} else {
  results[[length(results) + 1]] <-
    uafr_status_row("Optional live categorate smoke test", "SKIP",
                    "Set UAFR_STUDENT_LIVE=1 or run with --live to include live database queries.")
}

summary_table <- do.call(rbind, results)
uafr_note("", log_file = log_file)
uafr_note("Acceptance test summary:", log_file = log_file)
uafr_print_status(summary_table, log_file = log_file)

failures <- summary_table$Check[summary_table$Status == "FAIL"]
uafr_note("", log_file = log_file)
if (length(failures) > 0) {
  uafr_note("Acceptance test result: FAIL. Review failed item(s): ",
            paste(failures, collapse = ", "), log_file = log_file)
  stop("Student acceptance test failed. See ", log_file, call. = FALSE)
}

uafr_note("Acceptance test result: PASS.", log_file = log_file)
uafr_note("Log file: ", log_file, log_file = log_file)
