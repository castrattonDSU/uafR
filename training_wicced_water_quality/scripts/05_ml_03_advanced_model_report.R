message("WiCCED advanced ML step 3: advanced model review report")

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(cmd_file) > 0) {
  dirname(normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = FALSE))
} else {
  getwd()
}
manual_dir <- dirname(script_dir)

args <- commandArgs(trailingOnly = TRUE)
beginner_dir <- if (length(args) >= 1) args[[1]] else file.path(manual_dir, "results", "ml")
advanced_dir <- if (length(args) >= 2) args[[2]] else file.path(manual_dir, "results", "ml_advanced")
required_beginner <- c("ml_feature_table.csv", "ml_feature_dictionary.csv", "ml_model_card.md")
required_advanced <- c("advanced_ml_validation_folds.csv",
                       "advanced_ml_split_comparison.csv",
                       "advanced_ml_feature_set_metrics.csv",
                       "advanced_ml_permutation_importance.csv",
                       "advanced_ml_sensitivity_warnings.csv")
missing_beginner <- required_beginner[!file.exists(file.path(beginner_dir, required_beginner))]
missing_advanced <- required_advanced[!file.exists(file.path(advanced_dir, required_advanced))]
if (length(missing_beginner) > 0) {
  stop("Missing beginner ML outputs: ", paste(missing_beginner, collapse = ", "),
       ". Run scripts/04_ml_01 through scripts/04_ml_04 first.", call. = FALSE)
}
if (length(missing_advanced) > 0) {
  stop("Missing advanced ML outputs: ", paste(missing_advanced, collapse = ", "),
       ". Run scripts/05_ml_01 and scripts/05_ml_02 first.", call. = FALSE)
}

features <- utils::read.csv(file.path(beginner_dir, "ml_feature_table.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
dictionary <- utils::read.csv(file.path(beginner_dir, "ml_feature_dictionary.csv"),
                              stringsAsFactors = FALSE)
validation <- utils::read.csv(file.path(advanced_dir, "advanced_ml_validation_folds.csv"),
                              stringsAsFactors = FALSE)
comparison <- utils::read.csv(file.path(advanced_dir, "advanced_ml_split_comparison.csv"),
                              stringsAsFactors = FALSE)
feature_metrics <- utils::read.csv(file.path(advanced_dir, "advanced_ml_feature_set_metrics.csv"),
                                   stringsAsFactors = FALSE)
importance <- utils::read.csv(file.path(advanced_dir, "advanced_ml_permutation_importance.csv"),
                              stringsAsFactors = FALSE)
warnings <- utils::read.csv(file.path(advanced_dir, "advanced_ml_sensitivity_warnings.csv"),
                            stringsAsFactors = FALSE)

best_reg <- feature_metrics[is.finite(feature_metrics$regression_mae), , drop = FALSE]
best_reg <- if (nrow(best_reg) > 0) best_reg[which.min(best_reg$regression_mae), , drop = FALSE] else NULL
best_class <- feature_metrics[is.finite(feature_metrics$classification_balanced_accuracy), , drop = FALSE]
best_class <- if (nrow(best_class) > 0) best_class[which.max(best_class$classification_balanced_accuracy), , drop = FALSE] else NULL

top_reg_importance <- importance[importance$target == "salinity_ppt" &
                                   is.finite(importance$mean_importance), , drop = FALSE]
top_reg_importance <- top_reg_importance[order(top_reg_importance$mean_importance,
                                               decreasing = TRUE), , drop = FALSE]
top_class_importance <- importance[importance$target == "high_salinity" &
                                     is.finite(importance$mean_importance), , drop = FALSE]
top_class_importance <- top_class_importance[order(top_class_importance$mean_importance,
                                                   decreasing = TRUE), , drop = FALSE]

site_holdout_done <- any(validation$split_design == "leave_one_site_out")
warning_rows <- validation[nzchar(validation$warning), , drop = FALSE]
conductivity_used <- any(grepl("conductivity", feature_metrics$feature_set))

decision_gate <- data.frame(
  gate = c("beginner_ladder_complete",
           "advanced_validation_complete",
           "site_holdout_inspected",
           "feature_set_sensitivity_complete",
           "conductivity_decision_reviewed",
           "claim_ready_for_real_project"),
  status = c("ready",
             "ready",
             ifelse(site_holdout_done, "ready", "needs_work"),
             "ready",
             ifelse(conductivity_used, "requires_human_review", "not_applicable"),
             "not_ready_for_environmental_claims"),
  evidence_file = c("ml_model_card.md",
                    "advanced_ml_validation_folds.csv",
                    "advanced_ml_validation_folds.csv",
                    "advanced_ml_feature_set_metrics.csv",
                    "advanced_ml_sensitivity_warnings.csv",
                    "advanced_ml_model_review.md"),
  required_action = c(
    "Student can explain the beginner feature table, split, baselines, predictions, and model card before using advanced outputs.",
    "Student compares time-holdout and blocked validation rather than relying on one split.",
    "Student checks whether model behavior changes when an entire site is held out.",
    "Student compares feature sets on the same rows and split before changing the preferred model.",
    "Conductivity may be scientifically valid in some designs, but it must be justified because it is closely related to salinity.",
    "Replace synthetic teaching data with approved project data, rerun QA/QC, rerun beginner and advanced scripts, and update the model card."
  ),
  stringsAsFactors = FALSE
)

recommendations <- data.frame(
  recommendation = c("keep_beginner_model_as_reference",
                     "use_blocked_validation_before_final_claims",
                     "treat_conductivity_as_sensitivity_not_default",
                     "prefer_interpretable_feature_sets",
                     "document_skipped_overfit_models",
                     "do_not_claim_causation_from_importance"),
  reason = c(
    "The beginner model is the baseline that every advanced model must beat under the same evaluation design.",
    "Water-quality data are structured by site and date; random or single holdout results can be too optimistic.",
    "Conductivity can strongly track salinity and may turn the model into a near-duplicate sensor rather than a broader predictor.",
    "Small datasets do not support large feature sets or highly flexible models without strong validation.",
    "A skipped model is useful evidence that the dataset is too small for that feature set.",
    "Permutation importance shows what the fitted model used under a split, not the physical cause of the response."
  ),
  stringsAsFactors = FALSE
)

format_row <- function(row, columns) {
  if (is.null(row) || nrow(row) == 0) return("- No finite result available.")
  paste0("- ", paste(paste(columns, row[1, columns], sep = ": "), collapse = "; "))
}

review <- c(
  "# Advanced ML Model Review",
  "",
  "## Prerequisite",
  "",
  "This report is an extension, not a replacement for the beginner ML ladder. It should be interpreted only after the student can explain the feature table, leakage exclusions, baseline comparison, held-out predictions, and beginner model card.",
  "",
  "## What This Extension Adds",
  "",
  "- Time-holdout and leave-one-site-out validation.",
  "- Feature-set sensitivity, including conductivity as a flagged sensitivity feature.",
  "- Permutation importance on held-out rows.",
  "- A decision gate that separates training workflow results from real project claims.",
  "",
  "## Teaching Dataset Status",
  "",
  paste("- Rows in feature table:", nrow(features)),
  paste("- Candidate predictor rows in dictionary:", sum(dictionary$role %in% c("predictor", "excluded_predictor", "context_predictor_optional"))),
  "- Data status: synthetic teaching data only.",
  "",
  "## Best Teaching-Split Metrics",
  "",
  "Best regression feature set by held-out MAE:",
  format_row(best_reg, c("feature_set", "regression_mae", "regression_test_r2")),
  "",
  "Best classification feature set by held-out balanced accuracy:",
  format_row(best_class, c("feature_set", "classification_balanced_accuracy", "classification_sensitivity", "classification_specificity")),
  "",
  "## Top Permutation Importance Signals",
  "",
  "Salinity regression:",
  if (nrow(top_reg_importance) > 0) {
    paste0("- ", head(top_reg_importance$feature_set, 5), " / ",
           head(top_reg_importance$feature, 5), ": ",
           signif(head(top_reg_importance$mean_importance, 5), 4),
           " mean held-out MAE increase")
  } else {
    "- No finite regression importance rows."
  },
  "",
  "High-salinity classification:",
  if (nrow(top_class_importance) > 0) {
    paste0("- ", head(top_class_importance$feature_set, 5), " / ",
           head(top_class_importance$feature, 5), ": ",
           signif(head(top_class_importance$mean_importance, 5), 4),
           " mean balanced-accuracy drop")
  } else {
    "- No finite classification importance rows."
  },
  "",
  "## Validation Warnings",
  "",
  if (nrow(warning_rows) > 0) {
    paste0("- ", unique(warning_rows$warning[nzchar(warning_rows$warning)]))
  } else {
    "- No model-fit warnings were recorded in the validation output."
  },
  "",
  "## Required Interpretation",
  "",
  "- A feature-set improvement is not automatically a scientific improvement.",
  "- Conductivity sensitivity must be justified by the intended prediction setting.",
  "- Site-holdout performance is often more important than a single time split when the model is expected to transfer across sites.",
  "- Permutation importance is not causal evidence.",
  "- Synthetic teaching outputs do not support real WiCCED environmental claims.",
  "",
  "## Files to Review",
  "",
  "- `advanced_ml_validation_folds.csv`",
  "- `advanced_ml_split_comparison.csv`",
  "- `advanced_ml_feature_set_metrics.csv`",
  "- `advanced_ml_permutation_importance.csv`",
  "- `advanced_ml_sensitivity_warnings.csv`",
  "- `advanced_ml_decision_gate.csv`",
  "- `advanced_ml_final_recommendations.csv`"
)

utils::write.csv(decision_gate,
                 file.path(advanced_dir, "advanced_ml_decision_gate.csv"),
                 row.names = FALSE)
utils::write.csv(recommendations,
                 file.path(advanced_dir, "advanced_ml_final_recommendations.csv"),
                 row.names = FALSE)
writeLines(review, file.path(advanced_dir, "advanced_ml_model_review.md"))

message("Decision gates: ", nrow(decision_gate))
message("Wrote advanced model review outputs to: ", normalizePath(advanced_dir))
