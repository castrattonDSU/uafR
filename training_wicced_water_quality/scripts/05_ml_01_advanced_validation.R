message("WiCCED advanced ML step 1: blocked validation")

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

required_cols <- c("sample_id", "site_id", "sample_date", "split",
                   "salinity_ppt", "high_salinity")
missing_cols <- setdiff(required_cols, names(features))
if (length(missing_cols) > 0) {
  stop("Feature table is missing required columns: ",
       paste(missing_cols, collapse = ", "), call. = FALSE)
}

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
  precision <- if ((tp + fp) == 0) NA_real_ else tp / (tp + fp)
  c(accuracy = mean(predicted == actual, na.rm = TRUE),
    sensitivity = sensitivity,
    specificity = specificity,
    balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    precision = precision,
    true_positive = tp,
    true_negative = tn,
    false_positive = fp,
    false_negative = fn)
}
blank_row <- function(split_design, fold_id, target, model, feature_set,
                      train_rows, test_rows, features_used, warning = "") {
  data.frame(
    split_design = split_design,
    fold_id = fold_id,
    target = target,
    model = model,
    feature_set = feature_set,
    train_rows = train_rows,
    test_rows = test_rows,
    features_used = paste(features_used, collapse = "; "),
    rmse = NA_real_,
    mae = NA_real_,
    test_r2 = NA_real_,
    accuracy = NA_real_,
    sensitivity = NA_real_,
    specificity = NA_real_,
    balanced_accuracy = NA_real_,
    precision = NA_real_,
    true_positive = NA_real_,
    true_negative = NA_real_,
    false_positive = NA_real_,
    false_negative = NA_real_,
    warning = warning,
    stringsAsFactors = FALSE
  )
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
evaluate_fold <- function(train, test, split_design, fold_id, feature_set_name, requested_features) {
  used <- usable_features(train, test, requested_features)
  rows <- list()
  train_rows <- nrow(train)
  test_rows <- nrow(test)
  if (test_rows < 1 || train_rows < 4) {
    warn <- "Skipped because fold has too few training or testing rows."
    return(rbind(
      blank_row(split_design, fold_id, "salinity_ppt", "mean_only_baseline",
                feature_set_name, train_rows, test_rows, character(), warn),
      blank_row(split_design, fold_id, "high_salinity", "majority_class_baseline",
                feature_set_name, train_rows, test_rows, character(), warn)
    ))
  }
  if (length(used) >= train_rows - 2) {
    warn <- paste0("Skipped model fit because ", length(used),
                   " predictors are too many for ", train_rows,
                   " training rows. This is an overfitting risk.")
  } else {
    warn <- ""
  }

  reg_base <- rep(mean(train$salinity_ppt, na.rm = TRUE), test_rows)
  reg_base_row <- blank_row(split_design, fold_id, "salinity_ppt",
                            "mean_only_baseline", feature_set_name,
                            train_rows, test_rows, character())
  reg_base_row$rmse <- rmse(test$salinity_ppt, reg_base)
  reg_base_row$mae <- mae(test$salinity_ppt, reg_base)
  reg_base_row$test_r2 <- r2(test$salinity_ppt, reg_base)
  rows[[length(rows) + 1]] <- reg_base_row

  reg_row <- blank_row(split_design, fold_id, "salinity_ppt",
                       "linear_model", feature_set_name,
                       train_rows, test_rows, used, warn)
  if (length(used) > 0 && warn == "") {
    reg_formula <- stats::as.formula(paste("salinity_ppt ~", paste(used, collapse = " + ")))
    reg_model <- suppressWarnings(stats::lm(reg_formula, data = train))
    reg_pred <- suppressWarnings(as.numeric(stats::predict(reg_model, newdata = test)))
    if (all(is.finite(reg_pred))) {
      reg_row$rmse <- rmse(test$salinity_ppt, reg_pred)
      reg_row$mae <- mae(test$salinity_ppt, reg_pred)
      reg_row$test_r2 <- r2(test$salinity_ppt, reg_pred)
    } else {
      reg_row$warning <- "Regression prediction returned non-finite values."
    }
  }
  rows[[length(rows) + 1]] <- reg_row

  majority_class <- as.numeric(names(sort(table(train$high_salinity), decreasing = TRUE)[1]))
  class_base_pred <- rep(majority_class, test_rows)
  class_base <- blank_row(split_design, fold_id, "high_salinity",
                          "majority_class_baseline", feature_set_name,
                          train_rows, test_rows, character())
  cm <- class_metrics(test$high_salinity, class_base_pred)
  for (nm in names(cm)) class_base[[nm]] <- unname(cm[[nm]])
  rows[[length(rows) + 1]] <- class_base

  class_row <- blank_row(split_design, fold_id, "high_salinity",
                         "logistic_model", feature_set_name,
                         train_rows, test_rows, used, warn)
  if (length(unique(train$high_salinity)) < 2) {
    class_row$warning <- "Skipped logistic model because training fold has one class."
  } else if (length(used) > 0 && warn == "") {
    class_formula <- stats::as.formula(paste("high_salinity ~", paste(used, collapse = " + ")))
    class_model <- suppressWarnings(stats::glm(class_formula, data = train,
                                               family = stats::binomial()))
    class_prob <- suppressWarnings(as.numeric(stats::predict(class_model, newdata = test,
                                                             type = "response")))
    if (all(is.finite(class_prob))) {
      class_pred <- ifelse(class_prob >= 0.5, 1, 0)
      cm <- class_metrics(test$high_salinity, class_pred)
      for (nm in names(cm)) class_row[[nm]] <- unname(cm[[nm]])
    } else {
      class_row$warning <- "Logistic prediction returned non-finite values."
    }
  }
  rows[[length(rows) + 1]] <- class_row
  do.call(rbind, rows)
}

folds <- list(
  list(split_design = "time_holdout",
       fold_id = "train_early_test_late",
       train = which(features$split == "train"),
       test = which(features$split == "test"))
)
for (site in sort(unique(features$site_id))) {
  folds[[length(folds) + 1]] <- list(
    split_design = "leave_one_site_out",
    fold_id = paste0("heldout_site_", site),
    train = which(features$site_id != site),
    test = which(features$site_id == site)
  )
}

validation_rows <- list()
for (fold in folds) {
  for (feature_set_name in names(feature_sets)) {
    validation_rows[[length(validation_rows) + 1]] <- evaluate_fold(
      train = features[fold$train, , drop = FALSE],
      test = features[fold$test, , drop = FALSE],
      split_design = fold$split_design,
      fold_id = fold$fold_id,
      feature_set_name = feature_set_name,
      requested_features = feature_sets[[feature_set_name]]
    )
  }
}
validation_folds <- do.call(rbind, validation_rows)

mean_or_na <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}
key_cols <- c("split_design", "target", "model", "feature_set")
keys <- unique(validation_folds[key_cols])
comparison <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
  idx <- Reduce(`&`, Map(`==`, validation_folds[key_cols], keys[i, ]))
  group <- validation_folds[idx, , drop = FALSE]
  data.frame(
    keys[i, , drop = FALSE],
    folds_evaluated = nrow(group),
    folds_with_warnings = sum(nzchar(group$warning)),
    mean_rmse = mean_or_na(group$rmse),
    mean_mae = mean_or_na(group$mae),
    mean_test_r2 = mean_or_na(group$test_r2),
    mean_accuracy = mean_or_na(group$accuracy),
    mean_sensitivity = mean_or_na(group$sensitivity),
    mean_specificity = mean_or_na(group$specificity),
    mean_balanced_accuracy = mean_or_na(group$balanced_accuracy),
    stringsAsFactors = FALSE
  )
}))

notes <- c(
  "# Advanced Validation Notes",
  "",
  "This extension should be run only after the beginner ML ladder is understood and documented.",
  "",
  "The teaching dataset is intentionally small. The advanced outputs demonstrate validation design, leakage sensitivity, and model review rather than final environmental inference.",
  "",
  "Split designs included:",
  "- `time_holdout`: earlier dates train the model and later dates test it.",
  "- `leave_one_site_out`: each site is held out once to test whether the model transfers across sites.",
  "",
  "Feature sets included:",
  "- `core_water_quality`: the beginner approved features.",
  "- `core_plus_conductivity_sensitivity`: adds conductivity only as a documented sensitivity test.",
  "- `core_plus_context_sensitivity`: attempts to add one-hot context variables; small folds may be skipped if this creates too many predictors.",
  "",
  "Interpretation rule: an advanced model is not better because it is more complex. It is better only if it improves held-out performance, survives leakage review, and remains scientifically interpretable."
)

utils::write.csv(validation_folds,
                 file.path(out_dir, "advanced_ml_validation_folds.csv"),
                 row.names = FALSE)
utils::write.csv(comparison,
                 file.path(out_dir, "advanced_ml_split_comparison.csv"),
                 row.names = FALSE)
writeLines(notes, file.path(out_dir, "advanced_ml_validation_notes.md"))

message("Validation rows: ", nrow(validation_folds))
message("Wrote advanced validation outputs to: ", normalizePath(out_dir))
