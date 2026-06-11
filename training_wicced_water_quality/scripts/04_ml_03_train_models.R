message("WiCCED ML step 3: train beginner baseline models")

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
feature_path <- file.path(in_dir, "ml_feature_table.csv")
dictionary_path <- file.path(in_dir, "ml_feature_dictionary.csv")
if (!file.exists(feature_path)) stop("Missing ml_feature_table.csv. Run 04_ml_01_build_feature_table.R first.", call. = FALSE)
if (!file.exists(dictionary_path)) stop("Missing ml_feature_dictionary.csv. Run 04_ml_01_build_feature_table.R first.", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

features <- utils::read.csv(feature_path, stringsAsFactors = FALSE, check.names = FALSE)
features$sample_date <- as.Date(features$sample_date)
dictionary <- utils::read.csv(dictionary_path, stringsAsFactors = FALSE)
model_features <- dictionary$variable[dictionary$include_in_beginner_model == "yes"]
model_features <- model_features[vapply(features[model_features], is.numeric, logical(1))]

train <- features[features$split == "train", , drop = FALSE]
test <- features[features$split == "test", , drop = FALSE]
if (nrow(train) < 4 || nrow(test) < 1) stop("Not enough train/test rows for beginner ML models.", call. = FALSE)

rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))
mae <- function(actual, predicted) mean(abs(actual - predicted), na.rm = TRUE)
r2 <- function(actual, predicted) {
  denom <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  if (denom == 0) return(NA_real_)
  1 - sum((actual - predicted)^2, na.rm = TRUE) / denom
}

reg_formula <- stats::as.formula(paste("salinity_ppt ~", paste(model_features, collapse = " + ")))
reg_model <- stats::lm(reg_formula, data = train)
reg_pred <- as.numeric(stats::predict(reg_model, newdata = test))
reg_baseline <- rep(mean(train$salinity_ppt, na.rm = TRUE), nrow(test))

regression_metrics <- data.frame(
  task = "regression_salinity_ppt",
  model = c("mean_only_baseline", "linear_model_beginner_features"),
  rmse = c(rmse(test$salinity_ppt, reg_baseline), rmse(test$salinity_ppt, reg_pred)),
  mae = c(mae(test$salinity_ppt, reg_baseline), mae(test$salinity_ppt, reg_pred)),
  test_r2 = c(r2(test$salinity_ppt, reg_baseline), r2(test$salinity_ppt, reg_pred)),
  train_rows = nrow(train),
  test_rows = nrow(test),
  features = c("training mean only", paste(model_features, collapse = "; ")),
  stringsAsFactors = FALSE
)

class_formula <- stats::as.formula(paste("high_salinity ~", paste(model_features, collapse = " + ")))
class_model <- suppressWarnings(stats::glm(class_formula, data = train, family = stats::binomial()))
class_prob <- suppressWarnings(as.numeric(stats::predict(class_model, newdata = test, type = "response")))
class_pred <- ifelse(class_prob >= 0.5, 1, 0)
majority_class <- as.numeric(names(sort(table(train$high_salinity), decreasing = TRUE)[1]))
majority_pred <- rep(majority_class, nrow(test))

classification_metrics_for <- function(actual, predicted, model_name) {
  tp <- sum(predicted == 1 & actual == 1)
  tn <- sum(predicted == 0 & actual == 0)
  fp <- sum(predicted == 1 & actual == 0)
  fn <- sum(predicted == 0 & actual == 1)
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  data.frame(
    task = "classification_high_salinity",
    model = model_name,
    threshold = ifelse(model_name == "logistic_model_beginner_features", 0.5, NA),
    accuracy = mean(predicted == actual),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    precision = precision,
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn,
    train_rows = nrow(train),
    test_rows = nrow(test),
    stringsAsFactors = FALSE
  )
}

classification_metrics <- rbind(
  classification_metrics_for(test$high_salinity, majority_pred, "majority_class_baseline"),
  classification_metrics_for(test$high_salinity, class_pred, "logistic_model_beginner_features")
)

predictions <- data.frame(
  sample_id = test$sample_id,
  site_id = test$site_id,
  sample_date = as.character(test$sample_date),
  land_context = test$land_context,
  observed_salinity_ppt = test$salinity_ppt,
  predicted_salinity_mean_baseline = reg_baseline,
  predicted_salinity_linear_model = reg_pred,
  regression_residual = test$salinity_ppt - reg_pred,
  observed_high_salinity = test$high_salinity,
  predicted_high_salinity_majority_baseline = majority_pred,
  predicted_high_salinity_probability = class_prob,
  predicted_high_salinity_logistic = class_pred,
  stringsAsFactors = FALSE
)

coef_table <- data.frame(
  model = "linear_model_beginner_features",
  term = names(stats::coef(reg_model)),
  estimate = as.numeric(stats::coef(reg_model)),
  stringsAsFactors = FALSE
)
coef_table <- rbind(coef_table, data.frame(
  model = "logistic_model_beginner_features",
  term = names(stats::coef(class_model)),
  estimate = as.numeric(stats::coef(class_model)),
  stringsAsFactors = FALSE
))

summary_lines <- c(
  "# Beginner ML Model Comparison Summary",
  "",
  "## What Was Fit",
  "",
  "- Regression target: salinity_ppt.",
  "- Classification target: high_salinity.",
  "- Beginner features:",
  paste0("  - ", model_features),
  "",
  "## Required Interpretation",
  "",
  "These models are trained on synthetic teaching data. They demonstrate workflow, split discipline, baseline comparison, and model-card writing. They do not support environmental claims."
)

utils::write.csv(regression_metrics, file.path(out_dir, "ml_regression_metrics.csv"), row.names = FALSE)
utils::write.csv(classification_metrics, file.path(out_dir, "ml_classification_metrics.csv"), row.names = FALSE)
utils::write.csv(predictions, file.path(out_dir, "ml_predictions.csv"), row.names = FALSE)
utils::write.csv(coef_table, file.path(out_dir, "ml_model_coefficients.csv"), row.names = FALSE)
writeLines(summary_lines, file.path(out_dir, "ml_model_comparison_summary.md"))

message("Wrote ML step 3 outputs to: ", normalizePath(out_dir))
print(regression_metrics)
print(classification_metrics)
