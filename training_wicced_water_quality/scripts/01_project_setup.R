project_root <- getwd()
message("Creating WiCCED training project skeleton in: ", project_root)

dirs <- c(
  "data_raw",
  "data_processed",
  "results",
  "figures",
  "models",
  "logs",
  "scripts",
  "reports",
  "references"
)

for (dir in dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }
}

if (!file.exists("README.md")) {
  writeLines(c(
    "# WiCCED Water-Quality uafR/ML Training Project",
    "",
    "## Project question",
    "Replace this paragraph with the approved Week 1 research question.",
    "",
    "## Data rules",
    "Raw data stay in data_raw/. Processed data are written to data_processed/. Do not overwrite raw data.",
    "",
    "## Script order",
    "1. scripts/00_setup_check.R",
    "2. scripts/02_water_quality_qaqc.R",
    "3. scripts/03_uafr_chemical_screening_template.R",
    "4. scripts/04_ml_01_build_feature_table.R",
    "5. scripts/04_ml_02_explore_features.R",
    "6. scripts/04_ml_03_train_models.R",
    "7. scripts/04_ml_04_diagnostics_and_model_card.R",
    "",
    "## Advanced extension after beginner ML review",
    "8. scripts/05_ml_01_advanced_validation.R",
    "9. scripts/05_ml_02_permutation_importance_and_sensitivity.R",
    "10. scripts/05_ml_03_advanced_model_report.R",
    "",
    "## Current status",
    "Initial project skeleton created."
  ), "README.md")
}

log_file <- file.path("logs", "daily_research_log.md")
if (!file.exists(log_file)) {
  writeLines(c(
    "# Daily Research Log",
    "",
    "## Entry template",
    "- Date:",
    "- Goal:",
    "- Files opened:",
    "- Scripts or commands run:",
    "- Outputs created:",
    "- Errors and fixes:",
    "- Decisions:",
    "- Evidence consulted:",
    "- Interpretation limit:",
    "- Next action:"
  ), log_file)
}

message("Project skeleton ready.")
