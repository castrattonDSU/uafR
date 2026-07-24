message("WiCCED ML step 4: diagnostics and model card")

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

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[[1]] else file.path(manual_dir, "results", "ml")
out_dir <- if (length(args) >= 2) args[[2]] else in_dir
required_files <- c("ml_feature_table.csv", "ml_feature_dictionary.csv",
                    "ml_train_test_split_plan.csv", "ml_regression_metrics.csv",
                    "ml_classification_metrics.csv", "ml_predictions.csv",
                    "ml_model_coefficients.csv")
missing <- required_files[!file.exists(file.path(in_dir, required_files))]
if (length(missing) > 0) {
  stop("Missing ML outputs: ", paste(missing, collapse = ", "),
       ". Run steps 1-3 first.", call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

feature_table <- utils::read.csv(file.path(in_dir, "ml_feature_table.csv"),
                                 stringsAsFactors = FALSE, check.names = FALSE)
dictionary <- utils::read.csv(file.path(in_dir, "ml_feature_dictionary.csv"),
                              stringsAsFactors = FALSE)
split_plan <- utils::read.csv(file.path(in_dir, "ml_train_test_split_plan.csv"),
                              stringsAsFactors = FALSE)
regression_metrics <- utils::read.csv(file.path(in_dir, "ml_regression_metrics.csv"),
                                      stringsAsFactors = FALSE)
classification_metrics <- utils::read.csv(file.path(in_dir, "ml_classification_metrics.csv"),
                                          stringsAsFactors = FALSE)
predictions <- utils::read.csv(file.path(in_dir, "ml_predictions.csv"),
                               stringsAsFactors = FALSE)
coefficients <- utils::read.csv(file.path(in_dir, "ml_model_coefficients.csv"),
                                stringsAsFactors = FALSE)

pdf(file.path(out_dir, "ml_diagnostic_plots.pdf"), width = 8, height = 6)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(2, 2))
plot(predictions$observed_salinity_ppt,
     predictions$predicted_salinity_linear_model,
     xlab = "Observed salinity_ppt",
     ylab = "Predicted salinity_ppt",
     main = "Observed vs predicted")
abline(0, 1, col = "red")
plot(predictions$predicted_salinity_linear_model,
     predictions$regression_residual,
     xlab = "Predicted salinity_ppt",
     ylab = "Residual",
     main = "Residuals")
abline(h = 0, col = "red")
barplot(table(predictions$observed_high_salinity,
              predictions$predicted_high_salinity_logistic),
        beside = TRUE,
        legend.text = TRUE,
        main = "Observed vs predicted class",
        xlab = "Predicted high_salinity",
        ylab = "Rows")
plot(predictions$predicted_high_salinity_probability,
     predictions$observed_high_salinity,
     xlab = "Predicted high-salinity probability",
     ylab = "Observed class",
     main = "Classification probabilities")
dev.off()

best_reg <- regression_metrics[which.min(regression_metrics$mae), ]
best_class <- classification_metrics[which.max(classification_metrics$balanced_accuracy), ]
diagnostic_summary <- data.frame(
  diagnostic = c("best_regression_by_mae", "best_classification_by_balanced_accuracy",
                 "largest_absolute_regression_residual", "rows_in_test_set"),
  value = c(best_reg$model,
            best_class$model,
            predictions$sample_id[which.max(abs(predictions$regression_residual))],
            nrow(predictions)),
  note = c("Lower MAE is better for the continuous target.",
           "Balanced accuracy is used because class balance matters.",
           "Inspect this row before trusting model interpretation.",
           "Small test sets are useful for training but weak for final claims."),
  stringsAsFactors = FALSE
)

features_used <- dictionary$variable[dictionary$include_in_beginner_model == "yes"]
model_card <- c(
  "# WiCCED Beginner ML Model Card",
  "",
  "## Project Question",
  "Can basic water-quality measurements in the synthetic teaching dataset predict salinity and high-salinity class under a time-based split?",
  "",
  "## Data",
  paste("- Rows in feature table:", nrow(feature_table)),
  paste("- Training rows:", split_plan$train_rows[[1]]),
  paste("- Testing rows:", split_plan$test_rows[[1]]),
  "- Data status: synthetic teaching data only.",
  "",
  "## Targets",
  "- Regression target: salinity_ppt.",
  "- Classification target: high_salinity.",
  "",
  "## Features",
  paste0("- ", features_used),
  "",
  "## Split Design",
  paste("- Split type:", split_plan$split_type[[1]]),
  paste("- Cutoff date:", split_plan$cutoff_date[[1]]),
  paste("- Leakage rule:", split_plan$leakage_rule[[1]]),
  "",
  "## Regression Metrics",
  capture.output(print(regression_metrics, row.names = FALSE)),
  "",
  "## Classification Metrics",
  capture.output(print(classification_metrics, row.names = FALSE)),
  "",
  "## Coefficients",
  capture.output(print(coefficients, row.names = FALSE)),
  "",
  "## Diagnostic Notes",
  capture.output(print(diagnostic_summary, row.names = FALSE)),
  "",
  "## Limits",
  "- These models are for training only.",
  "- Results do not prove causation.",
  "- Results do not support regulatory, management, or site-safety decisions.",
  "- External validation with real project data is required before research interpretation.",
  "",
  "## Next Step",
  "Rerun the same ladder on an approved project dataset, then compare whether the same variables remain useful under a defensible split."
)

interpretation_guide <- c(
  "# Beginner ML Interpretation Guide",
  "",
  "Use these prompts after opening the metrics, predictions, coefficients, and diagnostic plots.",
  "",
  "1. Did the candidate regression model beat the mean-only baseline on MAE?",
  "2. Which test row had the largest absolute residual?",
  "3. Did the logistic model beat the majority-class baseline on balanced accuracy?",
  "4. Which mistake would matter more for the project: false positive or false negative?",
  "5. Which coefficient is largest, and why should that not be treated as causal proof?",
  "6. Which feature was excluded because it could leak salinity information?",
  "7. What new data would be needed before using this model for a real WiCCED claim?"
)

utils::write.csv(diagnostic_summary, file.path(out_dir, "ml_diagnostic_summary.csv"), row.names = FALSE)
writeLines(model_card, file.path(out_dir, "ml_model_card.md"))
writeLines(interpretation_guide, file.path(out_dir, "ml_interpretation_guide.md"))

message("Wrote ML step 4 outputs to: ", normalizePath(out_dir))
