#!/usr/bin/env Rscript

file_arg = grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path = if (length(file_arg) > 0) {
  sub("^--file=", "", file_arg[[1]])
} else {
  file.path(getwd(), "tools", "run_plant_chemistry_panel.R")
}
repo_root = normalizePath(file.path(dirname(script_path), ".."),
                          winslash = "/", mustWork = TRUE)
Sys.setenv(UAFR_DEV_ROOT = repo_root)
source(file.path(repo_root, "inst", "scripts",
                 "run_plant_chemistry_panel.R"), chdir = FALSE)
