message("WiCCED ML step 1: build a beginner feature table")

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
input <- if (length(args) >= 1) args[[1]] else file.path(manual_dir, "data", "wicced_teaching_water_quality.csv")
out_dir <- if (length(args) >= 2) args[[2]] else file.path(manual_dir, "results", "ml")
if (!file.exists(input)) {
  alt <- file.path("training_wicced_water_quality", "data", "wicced_teaching_water_quality.csv")
  if (file.exists(alt)) input <- alt
}
if (!file.exists(input)) stop("Could not find water-quality input CSV.", call. = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wq <- utils::read.csv(input, stringsAsFactors = FALSE)
required <- c("sample_id", "site_id", "sample_date", "land_context",
              "salinity_ppt", "conductivity_ms_cm", "temp_c", "do_mg_l",
              "ph", "turbidity_ntu", "nitrate_mg_l", "orthophosphate_mg_l",
              "chlorophyll_ug_l", "high_salinity")
missing <- setdiff(required, names(wq))
if (length(missing) > 0) {
  stop("Input is missing required ML teaching columns: ",
       paste(missing, collapse = ", "), call. = FALSE)
}

wq$sample_date <- as.Date(wq$sample_date)
if (any(is.na(wq$sample_date))) stop("sample_date could not be parsed as Date.", call. = FALSE)

numeric_cols <- c("salinity_ppt", "conductivity_ms_cm", "temp_c", "do_mg_l",
                  "ph", "turbidity_ntu", "nitrate_mg_l",
                  "orthophosphate_mg_l", "chlorophyll_ug_l", "high_salinity")
for (col in numeric_cols) {
  wq[[col]] <- suppressWarnings(as.numeric(wq[[col]]))
  if (any(is.na(wq[[col]]))) stop("Column ", col, " contains missing or non-numeric values.", call. = FALSE)
}

wq <- wq[order(wq$sample_date, wq$site_id), ]
context_matrix <- stats::model.matrix(~ land_context - 1, data = wq)
colnames(context_matrix) <- sub("^land_context", "context_", colnames(context_matrix))

feature_table <- cbind(
  wq[, c("sample_id", "site_id", "sample_date", "land_context",
         "salinity_ppt", "high_salinity", "conductivity_ms_cm", "temp_c",
         "do_mg_l", "ph", "turbidity_ntu", "nitrate_mg_l",
         "orthophosphate_mg_l", "chlorophyll_ug_l")],
  as.data.frame(context_matrix, check.names = FALSE)
)

cutoff <- stats::median(feature_table$sample_date)
feature_table$split <- ifelse(feature_table$sample_date <= cutoff, "train", "test")
feature_table$split_reason <- paste0("time split: train dates <= ", cutoff)

dictionary <- data.frame(
  variable = names(feature_table),
  role = "metadata",
  unit_or_encoding = "not_applicable",
  source = basename(input),
  transformation = "none",
  include_in_beginner_model = "no",
  leakage_note = "not a model feature",
  beginner_note = "Used for identification or documentation.",
  stringsAsFactors = FALSE
)

set_role <- function(vars, role, unit, include, leakage, note) {
  idx <- dictionary$variable %in% vars
  dictionary$role[idx] <<- role
  dictionary$unit_or_encoding[idx] <<- unit
  dictionary$include_in_beginner_model[idx] <<- include
  dictionary$leakage_note[idx] <<- leakage
  dictionary$beginner_note[idx] <<- note
}

set_role(c("salinity_ppt"), "regression_target", "ppt", "no",
         "target; never use as predictor for salinity models",
         "Continuous target for the beginner regression model.")
set_role(c("high_salinity"), "classification_target", "0/1", "no",
         "target; never use as predictor for classification models",
         "Binary target for the beginner classification model.")
set_role(c("temp_c", "do_mg_l", "ph", "turbidity_ntu", "nitrate_mg_l",
           "orthophosphate_mg_l", "chlorophyll_ug_l"),
         "predictor", "see variable name", "yes",
         "available before prediction in the teaching workflow",
         "Core beginner model feature.")
set_role(c("conductivity_ms_cm"), "excluded_predictor", "mS/cm", "no",
         "closely related to salinity; use only in a documented sensitivity model",
         "Excluded from the beginner salinity/high-salinity models to teach leakage caution.")
set_role(grep("^context_", dictionary$variable, value = TRUE),
         "context_predictor_optional", "one_hot_0_1", "optional",
         "can memorize site context in small data",
         "Use only after discussing whether context variables match the model claim.")
set_role(c("sample_id", "site_id", "sample_date", "land_context", "split", "split_reason"),
         "metadata", "not_applicable", "no",
         "identifier or split metadata; do not include in beginner model",
         "Useful for checking rows, groups, and split design.")

split_plan <- data.frame(
  split_type = "time_based",
  cutoff_date = as.character(cutoff),
  train_rows = sum(feature_table$split == "train"),
  test_rows = sum(feature_table$split == "test"),
  target_regression = "salinity_ppt",
  target_classification = "high_salinity",
  beginner_features = paste(dictionary$variable[dictionary$include_in_beginner_model == "yes"],
                            collapse = ", "),
  leakage_rule = "Do not use targets, identifiers, split columns, or conductivity in the beginner model.",
  stringsAsFactors = FALSE
)

utils::write.csv(feature_table, file.path(out_dir, "ml_feature_table.csv"), row.names = FALSE)
utils::write.csv(dictionary, file.path(out_dir, "ml_feature_dictionary.csv"), row.names = FALSE)
utils::write.csv(split_plan, file.path(out_dir, "ml_train_test_split_plan.csv"), row.names = FALSE)
utils::write.csv(feature_table[feature_table$split == "train", ], file.path(out_dir, "ml_training_rows.csv"), row.names = FALSE)
utils::write.csv(feature_table[feature_table$split == "test", ], file.path(out_dir, "ml_testing_rows.csv"), row.names = FALSE)

message("Feature table rows: ", nrow(feature_table))
message("Training rows: ", split_plan$train_rows)
message("Testing rows: ", split_plan$test_rows)
message("Wrote ML step 1 outputs to: ", normalizePath(out_dir))
