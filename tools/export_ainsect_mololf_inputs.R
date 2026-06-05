#!/usr/bin/env Rscript

script_path = sub("^--file=", "",
                  grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = FALSE)

load_uafr = function(repo_root) {
  if (requireNamespace("uafR", quietly = TRUE) &&
      exists("exportAiNsectMolOlfInputs", asNamespace("uafR"),
             inherits = FALSE)) {
    suppressPackageStartupMessages(library(uafR))
    return(invisible(TRUE))
  }
  if (file.exists(file.path(repo_root, "DESCRIPTION")) &&
      requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(repo_root, quiet = TRUE)
    return(invisible(TRUE))
  }
  stop("Could not load uafR with exportAiNsectMolOlfInputs(). Install the ",
       "current package or run this script from the uafR repository with ",
       "devtools installed.", call. = FALSE)
}

parse_args = function(args) {
  out = list()
  i = 1L
  while (i <= length(args)) {
    arg = args[[i]]
    if (!grepl("^--", arg)) {
      stop("Unexpected argument: ", arg, call. = FALSE)
    }
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

required_arg = function(args, name) {
  value = args[[name]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop("Missing required argument `--", gsub("_", "-", name), "`.",
         call. = FALSE)
  }
  value
}

usage = function() {
  cat(
    "Usage:\n",
    "Rscript tools/export_ainsect_mololf_inputs.R \\\n",
    "  --chem-id-csv path/to/identity.csv \\\n",
    "  --chem-quant-csv path/to/quant.csv \\\n",
    "  --out-dir path/to/mololf_export \\\n",
    "  --cache-dir path/to/pubchem_cache \\\n",
    "  --profile ms \\\n",
    "  --throttle 0.2\n",
    sep = ""
  )
}

main = function() {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) == 0 || any(args %in% c("--help", "-h"))) {
    usage()
    quit(status = ifelse(length(args) == 0, 1L, 0L))
  }
  parsed = parse_args(args)
  load_uafr(repo_root)

  result = exportAiNsectMolOlfInputs(
    chem_id_csv = required_arg(parsed, "chem_id_csv"),
    chem_quant_csv = required_arg(parsed, "chem_quant_csv"),
    out_dir = required_arg(parsed, "out_dir"),
    cache_dir = required_arg(parsed, "cache_dir"),
    profile = if (!is.null(parsed$profile)) parsed$profile else "ms",
    throttle = if (!is.null(parsed$throttle)) {
      as.numeric(parsed$throttle)
    } else {
      0.2
    },
    min_relative_abundance = if (!is.null(parsed$min_relative_abundance)) {
      as.numeric(parsed$min_relative_abundance)
    } else {
      0
    },
    source_label = if (!is.null(parsed$source_label)) {
      parsed$source_label
    } else {
      "uafR_pubchem_eo_gcms"
    },
    uafR_run_id = parsed$uafr_run_id,
    write_pubchem_audit = !isTRUE(parsed$no_pubchem_audit)
  )

  cat("aiNsect molecular-olfaction export complete.\n")
  cat("Compounds: ", result$paths$compounds, "\n", sep = "")
  cat("Abundance: ", result$paths$abundance, "\n", sep = "")
  cat("Unresolved: ", result$paths$unresolved, "\n", sep = "")
  cat("Summary: ", result$paths$summary, "\n", sep = "")
}

main()
