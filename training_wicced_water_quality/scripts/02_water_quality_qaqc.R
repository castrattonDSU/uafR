args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1) args[[1]] else file.path("data", "wicced_teaching_water_quality.csv")
if (!file.exists(input)) {
  alt <- file.path("training_wicced_water_quality", "data", "wicced_teaching_water_quality.csv")
  if (file.exists(alt)) {
    input <- alt
  }
}
if (!file.exists(input)) {
  stop("Could not find teaching water-quality CSV. Provide a path as the first argument.")
}

out_dir <- if (length(args) >= 2) args[[2]] else "results"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

wq <- utils::read.csv(input, stringsAsFactors = FALSE)
required <- c("sample_id", "site_id", "sample_date", "salinity_ppt",
              "conductivity_ms_cm", "temp_c", "do_mg_l", "ph",
              "turbidity_ntu", "nitrate_mg_l", "chlorophyll_ug_l")
missing_cols <- setdiff(required, names(wq))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

range_rules <- list(
  salinity_ppt = c(0, 45),
  conductivity_ms_cm = c(0, 80),
  temp_c = c(-2, 40),
  do_mg_l = c(0, 20),
  ph = c(0, 14),
  turbidity_ntu = c(0, 1000),
  nitrate_mg_l = c(0, 100),
  orthophosphate_mg_l = c(0, 50),
  chlorophyll_ug_l = c(0, 500)
)

flags <- data.frame()
for (col in intersect(names(range_rules), names(wq))) {
  lim <- range_rules[[col]]
  values <- suppressWarnings(as.numeric(wq[[col]]))
  bad <- which(is.na(values) | values < lim[1] | values > lim[2])
  if (length(bad) > 0) {
    flags <- rbind(flags, data.frame(
      sample_id = wq$sample_id[bad],
      variable = col,
      raw_value = as.character(wq[[col]][bad]),
      flag = ifelse(is.na(values[bad]), "missing_or_non_numeric", "outside_training_range"),
      rule = paste0(col, " must be between ", lim[1], " and ", lim[2]),
      stringsAsFactors = FALSE
    ))
  }
}

summary_rows <- lapply(names(wq), function(col) {
  values <- suppressWarnings(as.numeric(wq[[col]]))
  numeric_col <- !all(is.na(values))
  data.frame(
    variable = col,
    n = length(wq[[col]]),
    missing = sum(is.na(wq[[col]]) | wq[[col]] == ""),
    numeric = numeric_col,
    min = if (numeric_col) min(values, na.rm = TRUE) else NA_real_,
    max = if (numeric_col) max(values, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)

utils::write.csv(wq, file.path(out_dir, "water_quality_cleaned.csv"), row.names = FALSE)
utils::write.csv(flags, file.path(out_dir, "water_quality_qaqc_flags.csv"), row.names = FALSE)
utils::write.csv(summary, file.path(out_dir, "water_quality_variable_summary.csv"), row.names = FALSE)

message("Rows read: ", nrow(wq))
message("QA/QC flags: ", nrow(flags))
message("Wrote outputs to: ", normalizePath(out_dir))
