message("WiCCED ML step 2: explore features before modeling")

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
dictionary <- utils::read.csv(dictionary_path, stringsAsFactors = FALSE)
model_features <- dictionary$variable[dictionary$include_in_beginner_model == "yes"]
numeric_model_features <- model_features[vapply(features[model_features], is.numeric, logical(1))]

summarize_numeric <- function(var) {
  values <- features[[var]]
  data.frame(
    variable = var,
    n = length(values),
    missing = sum(is.na(values)),
    mean = mean(values, na.rm = TRUE),
    sd = stats::sd(values, na.rm = TRUE),
    min = min(values, na.rm = TRUE),
    median = stats::median(values, na.rm = TRUE),
    max = max(values, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}
summary_table <- do.call(rbind, lapply(c("salinity_ppt", "high_salinity", numeric_model_features), summarize_numeric))

cor_vars <- c("salinity_ppt", numeric_model_features)
correlation <- stats::cor(features[, cor_vars], use = "pairwise.complete.obs")
correlation_table <- data.frame(
  feature = rownames(correlation),
  correlation_with_salinity = as.numeric(correlation[, "salinity_ppt"]),
  stringsAsFactors = FALSE
)
correlation_table <- correlation_table[order(abs(correlation_table$correlation_with_salinity),
                                             decreasing = TRUE), ]

class_balance <- as.data.frame(table(features$split, features$high_salinity),
                               stringsAsFactors = FALSE)
names(class_balance) <- c("split", "high_salinity", "rows")

pdf(file.path(out_dir, "ml_exploration_plots.pdf"), width = 8, height = 6)
old_par <- par(no.readonly = TRUE)
on.exit(par(old_par), add = TRUE)
par(mfrow = c(2, 2))
hist(features$salinity_ppt, main = "Salinity target", xlab = "salinity_ppt")
boxplot(salinity_ppt ~ high_salinity, data = features,
        main = "Salinity by class", xlab = "high_salinity", ylab = "salinity_ppt")
plot(features$temp_c, features$salinity_ppt,
     xlab = "temp_c", ylab = "salinity_ppt", main = "Temperature vs salinity")
plot(features$do_mg_l, features$salinity_ppt,
     xlab = "do_mg_l", ylab = "salinity_ppt", main = "Dissolved oxygen vs salinity")
par(mfrow = c(2, 2))
plot(features$turbidity_ntu, features$salinity_ppt,
     xlab = "turbidity_ntu", ylab = "salinity_ppt", main = "Turbidity vs salinity")
plot(features$nitrate_mg_l, features$salinity_ppt,
     xlab = "nitrate_mg_l", ylab = "salinity_ppt", main = "Nitrate vs salinity")
plot(features$chlorophyll_ug_l, features$salinity_ppt,
     xlab = "chlorophyll_ug_l", ylab = "salinity_ppt", main = "Chlorophyll vs salinity")
barplot(tapply(features$high_salinity, features$split, sum),
        main = "High-salinity rows by split", ylab = "rows")
dev.off()

questions <- c(
  "# Beginner ML Exploration Questions",
  "",
  "1. Which feature has the strongest simple correlation with salinity?",
  "2. Does that feature make scientific sense, or could it be a proxy for something else?",
  "3. Are both high-salinity and lower-salinity rows present in training and testing?",
  "4. Which feature would be most dangerous to include because it may leak the target?",
  "5. Which plot shows the clearest pattern?",
  "6. Which plot makes the model look uncertain or limited?",
  "",
  "Required conclusion: exploration can suggest patterns, but it does not prove causation."
)

utils::write.csv(summary_table, file.path(out_dir, "ml_feature_summary.csv"), row.names = FALSE)
utils::write.csv(correlation_table, file.path(out_dir, "ml_correlation_with_salinity.csv"), row.names = FALSE)
utils::write.csv(as.data.frame(correlation), file.path(out_dir, "ml_correlation_matrix.csv"))
utils::write.csv(class_balance, file.path(out_dir, "ml_class_balance.csv"), row.names = FALSE)
writeLines(questions, file.path(out_dir, "ml_exploration_questions.md"))

message("Wrote ML step 2 outputs to: ", normalizePath(out_dir))
