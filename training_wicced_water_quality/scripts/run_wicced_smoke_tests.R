locate_manual_dir <- function() {
  cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  candidates <- character()
  if (length(cmd_file) > 0) {
    script_path <- normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = FALSE)
    candidates <- c(candidates, dirname(dirname(script_path)), dirname(script_path))
  }
  candidates <- c(
    candidates,
    getwd(),
    dirname(getwd()),
    file.path(getwd(), "training_wicced_water_quality")
  )
  for (candidate in unique(candidates)) {
    if (dir.exists(file.path(candidate, "scripts")) &&
        dir.exists(file.path(candidate, "data")) &&
        file.exists(file.path(candidate, "README.md"))) {
      return(normalizePath(candidate, mustWork = FALSE))
    }
  }
  normalizePath(getwd(), mustWork = FALSE)
}

manual_dir <- locate_manual_dir()
script_dir <- file.path(manual_dir, "scripts")

script_path <- function(name) {
  file.path(script_dir, name)
}

data_path <- function(name) {
  file.path(manual_dir, "data", name)
}

run_script <- function(path, args = character()) {
  if (is.na(path) || !file.exists(path)) stop("Missing script: ", path)
  old <- commandArgs
  message("Running ", path)
  system2(file.path(R.home("bin"), "Rscript"), c(path, args), stdout = TRUE, stderr = TRUE)
}

tmp <- tempfile("wicced_training_smoke_")
dir.create(tmp)
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(tmp)

dir.create("data", recursive = TRUE)
invisible(file.copy(
  data_path("wicced_teaching_water_quality.csv"),
  file.path("data", "wicced_teaching_water_quality.csv")
))

scripts <- c("00_setup_check.R", "01_project_setup.R", "02_water_quality_qaqc.R",
             "03_uafr_chemical_screening_template.R",
             "04_ml_01_build_feature_table.R",
             "04_ml_02_explore_features.R",
             "04_ml_03_train_models.R",
             "04_ml_04_diagnostics_and_model_card.R",
             "05_ml_01_advanced_validation.R",
             "05_ml_02_permutation_importance_and_sensitivity.R",
             "05_ml_03_advanced_model_report.R")
paths <- vapply(scripts, script_path, character(1))

output <- list()
for (path in paths[1:2]) {
  output[[basename(path)]] <- run_script(path)
}
output[["02_water_quality_qaqc.R"]] <- run_script(paths[[3]], c(file.path("data", "wicced_teaching_water_quality.csv"), "results"))
output[["03_uafr_chemical_screening_template.R"]] <- run_script(paths[[4]], c(
  data_path("wicced_simulated_gcms_hits.csv"),
  data_path("wicced_simulated_gcms_queries.csv"),
  data_path("wicced_simulated_gcms_compound_reference.csv"),
  data_path("wicced_simulated_gcms_sample_metadata.csv"),
  "results",
  "offline"
))
output[["04_ml_01_build_feature_table.R"]] <- run_script(paths[[5]], c(file.path("data", "wicced_teaching_water_quality.csv"), "models"))
output[["04_ml_02_explore_features.R"]] <- run_script(paths[[6]], c("models", "models"))
output[["04_ml_03_train_models.R"]] <- run_script(paths[[7]], c("models", "models"))
output[["04_ml_04_diagnostics_and_model_card.R"]] <- run_script(paths[[8]], c("models", "models"))
output[["05_ml_01_advanced_validation.R"]] <- run_script(paths[[9]], c("models", "models_advanced"))
output[["05_ml_02_permutation_importance_and_sensitivity.R"]] <- run_script(paths[[10]], c("models", "models_advanced"))
output[["05_ml_03_advanced_model_report.R"]] <- run_script(paths[[11]], c("models", "models_advanced"))

expected <- c(
  "setup_package_status.csv",
  "README.md",
  file.path("logs", "daily_research_log.md"),
  file.path("results", "water_quality_cleaned.csv"),
  file.path("results", "water_quality_qaqc_flags.csv"),
  file.path("results", "water_quality_variable_summary.csv"),
  file.path("results", "wicced_simulated_spread.rds"),
  file.path("results", "wicced_simulated_exact_matches.csv"),
  file.path("results", "wicced_simulated_long_abundance.csv"),
  file.path("results", "wicced_simulated_sample_group_summary.csv"),
  file.path("results", "wicced_simulated_compound_evidence_notes.csv"),
  file.path("results", "wicced_simulated_gcms_pipeline_summary.txt"),
  file.path("models", "ml_feature_table.csv"),
  file.path("models", "ml_feature_dictionary.csv"),
  file.path("models", "ml_train_test_split_plan.csv"),
  file.path("models", "ml_training_rows.csv"),
  file.path("models", "ml_testing_rows.csv"),
  file.path("models", "ml_feature_summary.csv"),
  file.path("models", "ml_correlation_with_salinity.csv"),
  file.path("models", "ml_correlation_matrix.csv"),
  file.path("models", "ml_class_balance.csv"),
  file.path("models", "ml_exploration_plots.pdf"),
  file.path("models", "ml_exploration_questions.md"),
  file.path("models", "ml_regression_metrics.csv"),
  file.path("models", "ml_classification_metrics.csv"),
  file.path("models", "ml_predictions.csv"),
  file.path("models", "ml_model_coefficients.csv"),
  file.path("models", "ml_model_comparison_summary.md"),
  file.path("models", "ml_diagnostic_plots.pdf"),
  file.path("models", "ml_diagnostic_summary.csv"),
  file.path("models", "ml_model_card.md"),
  file.path("models", "ml_interpretation_guide.md"),
  file.path("models_advanced", "advanced_ml_validation_folds.csv"),
  file.path("models_advanced", "advanced_ml_split_comparison.csv"),
  file.path("models_advanced", "advanced_ml_validation_notes.md"),
  file.path("models_advanced", "advanced_ml_feature_set_metrics.csv"),
  file.path("models_advanced", "advanced_ml_permutation_importance.csv"),
  file.path("models_advanced", "advanced_ml_sensitivity_warnings.csv"),
  file.path("models_advanced", "advanced_ml_sensitivity_plots.pdf"),
  file.path("models_advanced", "advanced_ml_decision_gate.csv"),
  file.path("models_advanced", "advanced_ml_final_recommendations.csv"),
  file.path("models_advanced", "advanced_ml_model_review.md")
)

missing <- expected[!file.exists(expected)]
if (length(missing) > 0) {
  stop("Smoke test missing expected outputs: ", paste(missing, collapse = ", "))
}

message("WiCCED training smoke test passed.")
message("Temporary project path: ", tmp)
