# Offline-friendly smoke tests for the training support scripts.

message("Sourcing training scripts")
source("training/scripts/00_install_check.R")
source("training/scripts/01_project_setup.R")
source("training/scripts/04_core_workflow.R")
source("training/scripts/06_trait_matrix_analysis.R")
source("training/scripts/07_visualization.R")
source("training/scripts/09_reproducibility_check.R")

message("Checking project setup in a temporary directory")
tmp_project <- file.path(tempdir(), paste0("uafR-training-", Sys.getpid()))
create_training_project(tmp_project)
project_report <- check_training_project(tmp_project)
stopifnot(all(project_report$Directories$Exists))
stopifnot(all(project_report$Files$Exists))

run_live <- identical(Sys.getenv("UAFR_TRAINING_LIVE"), "1")

if (requireNamespace("uafR", quietly = TRUE)) {
  message("Running offline core workflow with package example data")
  core <- run_core_workflow(live_lookup = FALSE)
  stopifnot(is.list(core))
  stopifnot("exact" %in% names(core))
} else {
  message("uafR is not installed; skipping package-dependent offline workflow smoke test")
}

if (run_live && requireNamespace("uafR", quietly = TRUE)) {
  message("Running service-dependent core workflow with package example data")
  core <- run_core_workflow(live_lookup = TRUE)
  stopifnot(is.list(core))
  stopifnot("exact" %in% names(core))
} else if (run_live) {
  message("uafR is not installed; skipping service-dependent workflow smoke test")
} else {
  message("Skipping service-dependent uafR workflow smoke test")
  message("Set UAFR_TRAINING_LIVE=1 to include it intentionally")
}

message("Training smoke tests completed")
