message("WiCCED advanced ML step 2: permutation importance and sensitivity")

cmd_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(cmd_file) > 0) {
  dirname(normalizePath(sub("^--file=", "", cmd_file[[1]]), mustWork = FALSE))
} else {
  getwd()
}
manual_dir <- dirname(script_dir)

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[[1]] else file.path(manual_dir, "results", "ml")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(manual_dir, "results", "ml_advanced")
required_files <- c("ml_feature_table.csv", "ml_feature_dictionary.csv",
                    "ml_train_test_split_plan.csv", "ml_model_card.md")
missing <- required_files[!file.exists(file.path(in_dir, required_files))]
if (length(missing) > 0) {
  stop("Missing beginner ML outputs: ", paste(missing, collapse = ", "),
       ". Run scripts/04_ml_01 through scripts/04_ml_04 before the advanced extension.",
       call. = FALSE)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

features <- utils::read.csv(file.path(in_dir, "ml_feature_table.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
dictionary <- utils::read.csv(file.path(in_dir, "ml_feature_dictionary.csv"),
                              stringsAsFactors = FALSE)
features$sample_date <- as.Date(features$sample_date)

core_features <- dictionary$variable[dictionary$include_in_beginner_model == "yes"]
core_features <- core_features[core_features %in% names(features)]
core_features <- core_features[vapply(features[core_features], is.numeric, logical(1))]
if (length(core_features) == 0) stop("No numeric beginner features found.", call. = FALSE)

feature_sets <- list(core_water_quality = core_features)
if ("conductivity_ms_cm" %in% names(features)) {
  feature_sets$core_plus_conductivity_sensitivity <- unique(c(core_features, "conductivity_ms_cm"))
}
context_features <- grep("^context_", names(features), value = TRUE)
if (length(context_features) > 0) {
  feature_sets$core_plus_context_sensitivity <- unique(c(core_features, context_features))
}

train <- features[features$split == "train", , drop = FALSE]
test <- features[features$split == "test", , drop = FALSE]
if (nrow(train) < 4 || nrow(test) < 1) {
  stop("Not enough train/test rows for advanced sensitivity checks.", call. = FALSE)
}

rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm = TRUE))
mae <- function(actual, predicted) mean(abs(actual - predicted), na.rm = TRUE)
r2 <- function(actual, predicted) {
  denom <- sum((actual - mean(actual, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  1 - sum((actual - predicted)^2, na.rm = TRUE) / denom
}
class_metrics <- function(actual, predicted) {
  tp <- sum(predicted == 1 & actual == 1, na.rm = TRUE)
  tn <- sum(predicted == 0 & actual == 0, na.rm = TRUE)
  fp <- sum(predicted == 1 & actual == 0, na.rm = TRUE)
  fn <- sum(predicted == 0 & actual == 1, na.rm = TRUE)
  sensitivity <- if ((tp + fn) == 0) NA_real_ else tp / (tp + fn)
  specificity <- if ((tn + fp) == 0) NA_real_ else tn / (tn + fp)
  c(accuracy = mean(predicted == actual, na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn)
}
usable_features <- function(train, test, requested) {
  requested <- requested[requested %in% names(train)]
  requested <- requested[vapply(train[requested], is.numeric, logical(1))]
  requested[vapply(requested, function(col) {
    all(is.finite(train[[col]])) &&
      all(is.finite(test[[col]])) &&
      stats::sd(train[[col]], na.rm = TRUE) > 0
  }, logical(1))]
}
fit_feature_set <- function(feature_set_name, requested_features) {
  used <- usable_features(train, test, requested_features)
  warning <- ""
  if (length(used) >= nrow(train) - 2) {
    warning <- paste0("Skipped because ", length(used),
                      " predictors are too many for ", nrow(train),
                      " training rows in the teaching split.")
  }
  out <- data.frame(
    feature_set = feature_set_name,
    train_rows = nrow(train),
    test_rows = nrow(test),
    features_used = paste(used, collapse = "; "),
    regression_rmse = NA_real_,
    regression_mae = NA_real_,
    regression_test_r2 = NA_real_,
    classification_accuracy = NA_real_,
    classification_sensitivity = NA_real_,
    classification_specificity = NA_real_,
    classification_balanced_accuracy = NA_real_,
    false_positive = NA_real_,
    false_negative = NA_real_,
    warning = warning,
    stringsAsFactors = FALSE
  )
  if (length(used) == 0 || warning != "") return(out)

  reg_formula <- stats::as.formula(paste("salinity_ppt ~", paste(used, collapse = " + ")))
  reg_model <- suppressWarnings(stats::lm(reg_formula, data = train))
  reg_pred <- suppressWarnings(as.numeric(stats::predict(reg_model, newdata = test)))
  if (all(is.finite(reg_pred))) {
    out$regression_rmse <- rmse(test$salinity_ppt, reg_pred)
    out$regression_mae <- mae(test$salinity_ppt, reg_pred)
    out$regression_test_r2 <- r2(test$salinity_ppt, reg_pred)
  } else {
    out$warning <- paste(out$warning, "Regression prediction returned non-finite values.")
  }

  if (length(unique(train$high_salinity)) >= 2) {
    class_formula <- stats::as.formula(paste("high_salinity ~", paste(used, collapse = " + ")))
    class_model <- suppressWarnings(stats::glm(class_formula, data = train,
                                               family = stats::binomial()))
    class_prob <- suppressWarnings(as.numeric(stats::predict(class_model, newdata = test,
                                                             type = "response")))
    if (all(is.finite(class_prob))) {
      class_pred <- ifelse(class_prob >= 0.5, 1, 0)
      cm <- class_metrics(test$high_salinity, class_pred)
      out$classification_accuracy <- unname(cm[["accuracy"]])
      out$classification_sensitivity <- unname(cm[["sensitivity"]])
      out$classification_specificity <- unname(cm[["specificity"]])
      out$classification_balanced_accuracy <- unname(cm[["balanced_accuracy"]])
      out$false_positive <- unname(cm[["false_positive"]])
      out$false_negative <- unname(cm[["false_negative"]])
    } else {
      out$warning <- paste(out$warning, "Logistic prediction returned non-finite values.")
    }
  } else {
    out$warning <- paste(out$warning, "Training split contains one class.")
  }
  out
}

feature_set_metrics <- do.call(rbind, Map(fit_feature_set, names(feature_sets), feature_sets))

set.seed(20260604)
permutation_importance <- function(feature_set_name, requested_features, reps = 50) {
  used <- usable_features(train, test, requested_features)
  if (length(used) == 0 || length(used) >= nrow(train) - 2) {
    return(data.frame(
      feature_set = feature_set_name,
      target = character(),
      feature = character(),
      importance_metric = character(),
      mean_importance = numeric(),
      sd_importance = numeric(),
      repetitions = integer(),
      stringsAsFactors = FALSE
    ))
  }
  reg_formula <- stats::as.formula(paste("salinity_ppt ~", paste(used, collapse = " + ")))
  reg_model <- suppressWarnings(stats::lm(reg_formula, data = train))
  reg_base_pred <- suppressWarnings(as.numeric(stats::predict(reg_model, newdata = test)))
  reg_base_mae <- mae(test$salinity_ppt, reg_base_pred)

  class_model <- NULL
  class_base_ba <- NA_real_
  if (length(unique(train$high_salinity)) >= 2) {
    class_formula <- stats::as.formula(paste("high_salinity ~", paste(used, collapse = " + ")))
    class_model <- suppressWarnings(stats::glm(class_formula, data = train,
                                               family = stats::binomial()))
    class_prob <- suppressWarnings(as.numeric(stats::predict(class_model, newdata = test,
                                                             type = "response")))
    if (all(is.finite(class_prob))) {
      class_base_ba <- unname(class_metrics(test$high_salinity,
                                            ifelse(class_prob >= 0.5, 1, 0))[["balanced_accuracy"]])
    }
  }

  rows <- list()
  for (feature in used) {
    reg_changes <- numeric(reps)
    class_changes <- numeric(reps)
    for (i in seq_len(reps)) {
      perm_test <- test
      perm_test[[feature]] <- sample(perm_test[[feature]], replace = FALSE)
      perm_reg_pred <- suppressWarnings(as.numeric(stats::predict(reg_model, newdata = perm_test)))
      reg_changes[[i]] <- mae(test$salinity_ppt, perm_reg_pred) - reg_base_mae
      if (!is.null(class_model) && is.finite(class_base_ba)) {
        perm_prob <- suppressWarnings(as.numeric(stats::predict(class_model, newdata = perm_test,
                                                                type = "response")))
        if (all(is.finite(perm_prob))) {
          perm_ba <- unname(class_metrics(test$high_salinity,
                                          ifelse(perm_prob >= 0.5, 1, 0))[["balanced_accuracy"]])
          class_changes[[i]] <- class_base_ba - perm_ba
        } else {
          class_changes[[i]] <- NA_real_
        }
      } else {
        class_changes[[i]] <- NA_real_
      }
    }
    rows[[length(rows) + 1]] <- data.frame(
      feature_set = feature_set_name,
      target = "salinity_ppt",
      feature = feature,
      importance_metric = "heldout_mae_increase_when_permuted",
      mean_importance = mean(reg_changes, na.rm = TRUE),
      sd_importance = stats::sd(reg_changes, na.rm = TRUE),
      repetitions = reps,
      stringsAsFactors = FALSE
    )
    rows[[length(rows) + 1]] <- data.frame(
      feature_set = feature_set_name,
      target = "high_salinity",
      feature = feature,
      importance_metric = "heldout_balanced_accuracy_drop_when_permuted",
      mean_importance = mean(class_changes, na.rm = TRUE),
      sd_importance = stats::sd(class_changes, na.rm = TRUE),
      repetitions = reps,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

importance_sets <- feature_sets[names(feature_sets) %in%
                                  c("core_water_quality",
                                    "core_plus_conductivity_sensitivity")]
importance <- do.call(rbind, Map(permutation_importance, names(importance_sets), importance_sets))

warnings <- data.frame(
  issue = c("advanced_extension_prerequisite",
            "conductivity_sensitivity",
            "context_feature_sensitivity",
            "small_teaching_dataset",
            "permutation_importance_limit"),
  guidance = c(
    "Run and understand the four beginner ML scripts before interpreting these outputs.",
    "If conductivity improves performance, that may reflect its physical relationship to salinity rather than a generalizable independent predictor.",
    "Context features can memorize site type in small datasets; use them only when the project claim allows that context at prediction time.",
    "The synthetic teaching dataset is small. Use these outputs to learn the workflow, not to make environmental claims.",
    "Permutation importance shows model dependence under the chosen split. It is not causal evidence."
  ),
  stringsAsFactors = FALSE
)

pdf(file.path(out_dir, "advanced_ml_sensitivity_plots.pdf"), width = 8, height = 6)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(2, 2))
barplot(feature_set_metrics$regression_mae,
        names.arg = feature_set_metrics$feature_set,
        las = 2,
        ylab = "Held-out MAE",
        main = "Feature-set regression sensitivity")
barplot(feature_set_metrics$classification_balanced_accuracy,
        names.arg = feature_set_metrics$feature_set,
        las = 2,
        ylab = "Balanced accuracy",
        main = "Feature-set classification sensitivity")
reg_imp <- importance[importance$target == "salinity_ppt", , drop = FALSE]
reg_imp <- reg_imp[order(reg_imp$mean_importance, decreasing = TRUE), , drop = FALSE]
if (nrow(reg_imp) > 0) {
  barplot(head(reg_imp$mean_importance, 8),
          names.arg = head(reg_imp$feature, 8),
          las = 2,
          ylab = "MAE increase",
          main = "Permutation importance: salinity")
} else {
  plot.new()
  title("No regression importance rows")
}
class_imp <- importance[importance$target == "high_salinity", , drop = FALSE]
class_imp <- class_imp[order(class_imp$mean_importance, decreasing = TRUE), , drop = FALSE]
if (nrow(class_imp) > 0) {
  barplot(head(class_imp$mean_importance, 8),
          names.arg = head(class_imp$feature, 8),
          las = 2,
          ylab = "Balanced accuracy drop",
          main = "Permutation importance: class")
} else {
  plot.new()
  title("No classification importance rows")
}
dev.off()

utils::write.csv(feature_set_metrics,
                 file.path(out_dir, "advanced_ml_feature_set_metrics.csv"),
                 row.names = FALSE)
utils::write.csv(importance,
                 file.path(out_dir, "advanced_ml_permutation_importance.csv"),
                 row.names = FALSE)
utils::write.csv(warnings,
                 file.path(out_dir, "advanced_ml_sensitivity_warnings.csv"),
                 row.names = FALSE)

message("Feature sets evaluated: ", nrow(feature_set_metrics))
message("Permutation importance rows: ", nrow(importance))
message("Wrote advanced sensitivity outputs to: ", normalizePath(out_dir))
