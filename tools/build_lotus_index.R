#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "build_lotus_index.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

load_uafr = function(repo_root) {
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  if (requireNamespace("uafR", quietly = TRUE) &&
      exists("buildLotusIndex", asNamespace("uafR"), inherits = FALSE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  stop("Could not load uafR with buildLotusIndex(). Install the current ",
       "package or run this script from the uafR repository with devtools ",
       "installed.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) stop("Unexpected argument: ", arg, call. = FALSE)
    key = sub("^--", "", arg)
    value = TRUE
    if (grepl("=", key, fixed = TRUE)) {
      parts = strsplit(key, "=", fixed = TRUE)[[1]]
      key = parts[[1]]
      value = paste(parts[-1], collapse = "=")
    } else if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
      value = args[[i + 1L]]
      i = i + 1L
    }
    key = gsub("-", "_", key)
    out[[key]] = value
    i = i + 1L
  }
  out
}

usage = function() {
  cat(
    "Usage:\n",
    "Rscript tools/build_lotus_index.R \\\n",
    "  --input path/to/lotus_flat_export.jsonl \\\n",
    "  --out-file path/to/lotus_compact_index.csv \\\n",
    "  --overwrite\n\n",
    "Directory inputs are supported:\n",
    "Rscript tools/build_lotus_index.R \\\n",
    "  --input path/to/lotus_flat_exports/ \\\n",
    "  --out-file path/to/lotus_compact_index.rds \\\n",
    "  --format rds \\\n",
    "  --overwrite\n\n",
    "Notes:\n",
    "- Supported flat inputs: CSV, TSV, JSON, JSONL, NDJSON, and RDS.\n",
    "- Raw LOTUS MongoDB zip/BSON dumps can be flattened with ",
    "tools/flatten_lotus_mongo_dump.py. For production plant panels, use ",
    "that helper with --lookup-dir so uafR can query species/genus/family ",
    "shards without reading the complete compact CSV.\n",
    "- The output index can be passed to resolvePlantPhytochemistry(..., ",
    "lotus_index = 'lotus_compact_index.csv') or to the pilot wrapper with ",
    "--lotus-index.\n",
    sep = ""
  )
}

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop("Missing required argument `--", gsub("_", "-", name), "`.",
         call. = FALSE)
  }
  value
}

optional_arg = function(args, name, default = NULL) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    return(default)
  }
  value
}

flag_arg = function(args, name, default = FALSE) {
  value = args[[name]]
  if (is.null(value)) return(default)
  if (identical(value, TRUE)) return(TRUE)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

split_arg = function(value) {
  pieces = unlist(strsplit(value, "\\s*[,;]\\s*", perl = TRUE),
                  use.names = FALSE)
  pieces = trimws(pieces)
  pieces[nzchar(pieces)]
}

main = function() {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1L, 0L))
  }
  parsed = parse_args(args)
  load_uafr(repo_root)

  input = split_arg(required_arg(parsed, "input"))
  out_file = required_arg(parsed, "out_file")
  result = buildLotusIndex(
    input = input,
    out_file = out_file,
    format = optional_arg(parsed, "format", "auto"),
    overwrite = flag_arg(parsed, "overwrite", FALSE),
    manifest_file = optional_arg(parsed, "manifest_file", NULL),
    recursive = !flag_arg(parsed, "no_recursive", FALSE),
    strict = flag_arg(parsed, "strict", FALSE),
    progress = flag_arg(parsed, "progress", TRUE)
  )

  summary = result$BuildSummary
  cat("LOTUS compact index build complete.\n")
  cat("Output index: ", summary$out_file[[1]], "\n", sep = "")
  cat("Manifest: ", summary$manifest_file[[1]], "\n", sep = "")
  cat("Input sources: ", summary$source_count[[1]], "\n", sep = "")
  cat("Index rows: ", summary$index_row_count[[1]], "\n", sep = "")
  cat("Species: ", summary$species_count[[1]], "\n", sep = "")
  cat("Compounds: ", summary$compound_count[[1]], "\n", sep = "")
  cat("Rows with references: ", summary$rows_with_reference[[1]], "\n",
      sep = "")
}

main()
