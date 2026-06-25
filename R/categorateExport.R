#' Export enriched categorate results for review and sharing
#'
#' @description
#' `exportCategorateWorkbook()` writes the most useful `categorate()` result
#' tables to a researcher-friendly export bundle. It always supports a directory
#' of CSV files with an `ExportManifest` table. If `openxlsx` or `writexl` is
#' installed, it can also write a multi-sheet `.xlsx` workbook.
#'
#' @param x A list returned by `categorate()`, preferably with
#' `detail = "research"` or `detail = "full"`.
#' @param path Output path. For `format = "csv"`, this is a directory. For
#' `format = "xlsx"`, this is an `.xlsx` file.
#' @param tables Optional character vector of result table names to export. If
#' `NULL`, a curated analysis-ready set is exported.
#' @param format Export format: `"csv"`, `"xlsx"`, or `"auto"`. `"auto"` writes
#' `.xlsx` when `path` ends in `.xlsx` and an Excel writer is installed;
#' otherwise it writes a CSV bundle.
#' @param include_raw Logical. If `TRUE` and `tables = NULL`, export every data
#' frame in `x`, including raw source tables. If `FALSE`, export the curated
#' analysis and diagnostics tables.
#' @param include_empty Logical. If `TRUE`, include empty data frames so the
#' bundle preserves expected schema. If `FALSE`, omit empty tables.
#' @param overwrite Logical. If `TRUE`, replace an existing output file or
#' directory.
#' @param max_cell_chars Maximum characters retained in a single cell. Longer
#' values are truncated to keep spreadsheets responsive.
#'
#' @return A manifest data frame describing exported tables, row/column counts,
#' sheet names, file names, and output path.
#'
#' @examples
#' \dontrun{
#' result = categorate(compounds, library_data, detail = "research")
#' manifest = exportCategorateWorkbook(result, "categorate_export",
#'                                     format = "csv", overwrite = TRUE)
#' manifest
#' }
#'
#' @export
exportCategorateWorkbook = function(x,
                                    path,
                                    tables = NULL,
                                    format = c("auto", "xlsx", "csv"),
                                    include_raw = FALSE,
                                    include_empty = TRUE,
                                    overwrite = FALSE,
                                    max_cell_chars = 30000) {
  format = match.arg(format)
  if (missing(path) || is.null(path) || length(.uaf_non_empty(path)) < 1) {
    stop("`path` is required.", call. = FALSE)
  }
  if (!is.list(x) || is.data.frame(x)) {
    stop("`x` must be a categorate result list or a named list of data frames.",
         call. = FALSE)
  }

  resolved = .categorate_export_resolve_path(path, format)
  export_tables = .categorate_export_tables(
    x = x,
    tables = tables,
    include_raw = include_raw,
    include_empty = include_empty,
    max_cell_chars = max_cell_chars
  )
  if (length(export_tables) < 1) {
    stop("No exportable data frames were found.", call. = FALSE)
  }

  all_table_names = c("ExportManifest", names(export_tables))
  sheet_names = .categorate_export_sheet_names(all_table_names)
  file_names = .categorate_export_file_names(all_table_names)
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

  manifest = data.frame(
    Table = all_table_names,
    SheetName = sheet_names,
    FileName = if (resolved$format == "csv") file_names else NA_character_,
    RowCount = c(NA_integer_,
                 vapply(export_tables, nrow, integer(1))),
    ColumnCount = c(NA_integer_,
                    vapply(export_tables, ncol, integer(1))),
    Format = resolved$format,
    OutputPath = resolved$path,
    CreatedAt = created_at,
    stringsAsFactors = FALSE
  )
  manifest$RowCount[manifest$Table == "ExportManifest"] = nrow(manifest)
  manifest$ColumnCount[manifest$Table == "ExportManifest"] = ncol(manifest)

  full_tables = c(list(ExportManifest = manifest), export_tables)
  names(full_tables) = sheet_names

  if (resolved$format == "xlsx") {
    .categorate_write_xlsx(full_tables, resolved$path, overwrite)
  } else {
    .categorate_write_csv_bundle(full_tables, file_names, resolved$path,
                                 overwrite)
  }

  row.names(manifest) = NULL
  manifest
}

.categorate_export_tables = function(x, tables, include_raw, include_empty,
                                     max_cell_chars) {
  available = names(x)[vapply(x, is.data.frame, logical(1))]
  requested = .uaf_non_empty(tables)
  if (length(requested) > 0) {
    missing = setdiff(requested, available)
    if (length(missing) > 0) {
      warning("Requested table(s) not found and skipped: ",
              paste(missing, collapse = ", "), call. = FALSE)
    }
    table_names = intersect(requested, available)
  } else if (isTRUE(include_raw)) {
    table_names = available
  } else {
    table_names = intersect(.categorate_export_default_tables(), available)
  }

  out = list()
  for (table_name in table_names) {
    table = x[[table_name]]
    if (!isTRUE(include_empty) && nrow(table) < 1) next
    out[[table_name]] = .categorate_export_prepare_table(table,
                                                         max_cell_chars)
  }
  out
}

.categorate_export_default_tables = function() {
  c("ChemicalTraitReport", "DerivedGroups", "ChemicalTraitMatrix",
    "ChemicalTraitOntologyMatrix", "ChemicalTraitSummary",
    "ChemicalTraitSimilarity", "ChemicalTraits", "ChemicalTraitOntology",
    "ChemicalTraitEvidence", "ChemicalMeasurements",
    "ChemicalMeasurementSummary", "ChemicalHazards", "ChemicalUses",
    "ChemicalBioassays", "ChemicalBioactivities", "ChemicalTargets",
    "ChemicalPotencies", "ChemicalTaxonomy", "ChemicalOccurrences",
    "ChemicalPathwayRoles",
    "KEGGReactionParticipants", "SourceCoverage", "TableQuality",
    "SourceDiagnostics", "ValidationSummary", "ValidationIssues",
    "DataDictionary", "Provenance")
}

.categorate_export_prepare_table = function(table, max_cell_chars) {
  table = as.data.frame(table, stringsAsFactors = FALSE)
  for (col in colnames(table)) {
    if (is.list(table[[col]])) {
      table[[col]] = vapply(table[[col]], function(value) {
        .pubchem_collapse(unlist(value, use.names = FALSE))
      }, character(1))
    }
    if (is.factor(table[[col]])) {
      table[[col]] = as.character(table[[col]])
    }
    if (inherits(table[[col]], "POSIXt")) {
      table[[col]] = format(table[[col]], "%Y-%m-%dT%H:%M:%S%z")
    }
    if (is.character(table[[col]])) {
      table[[col]] = .categorate_export_clean_character(table[[col]])
      table[[col]] = .categorate_export_truncate(table[[col]],
                                                 max_cell_chars)
    }
  }
  row.names(table) = NULL
  table
}

.categorate_export_clean_character = function(x) {
  x = enc2utf8(as.character(x))
  keep = !is.na(x)
  if (any(keep)) {
    x[keep] = gsub("\r\n|\r|\n", " | ", x[keep], perl = TRUE)
    x[keep] = gsub("[\001-\010\013\014\016-\037\177]", " ",
                   x[keep], perl = TRUE)
  }
  x
}

.categorate_export_truncate = function(x, max_cell_chars) {
  max_cell_chars = suppressWarnings(as.integer(max_cell_chars[[1]]))
  if (is.na(max_cell_chars) || max_cell_chars < 100) max_cell_chars = 100
  too_long = !is.na(x) & nchar(x, type = "chars") > max_cell_chars
  x[too_long] = paste0(substr(x[too_long], 1, max_cell_chars), " ...")
  x
}

.categorate_export_resolve_path = function(path, format) {
  path = .uaf_non_empty(path)[[1]]
  has_xlsx_ext = grepl("\\.xlsx$", path, ignore.case = TRUE)
  xlsx_available = .categorate_xlsx_writer_available()
  if (format == "auto") {
    if (has_xlsx_ext && xlsx_available) {
      return(list(format = "xlsx", path = path))
    }
    if (has_xlsx_ext && !xlsx_available) {
      warning("No XLSX writer is installed. Install `openxlsx` or `writexl` ",
              "to write .xlsx files. Writing a CSV bundle instead.",
              call. = FALSE)
      path = sub("\\.xlsx$", "_csv", path, ignore.case = TRUE)
    }
    return(list(format = "csv", path = path))
  }
  if (format == "xlsx" && !xlsx_available) {
    stop("Writing .xlsx files requires the optional package `openxlsx` or ",
         "`writexl`. Install one of those packages or use `format = \"csv\"`.",
         call. = FALSE)
  }
  list(format = format, path = path)
}

.categorate_xlsx_writer_available = function() {
  any(vapply(.categorate_xlsx_writers(), function(package) {
    requireNamespace(package, quietly = TRUE)
  }, logical(1)))
}

.categorate_xlsx_writers = function() {
  c("openxlsx", "writexl")
}

.categorate_write_xlsx = function(tables, path, overwrite) {
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Output file already exists. Use `overwrite = TRUE` to replace it.",
         call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writers = .categorate_xlsx_writers()
  if (requireNamespace(writers[[1]], quietly = TRUE)) {
    create_workbook = getExportedValue(writers[[1]], "createWorkbook")
    add_worksheet = getExportedValue(writers[[1]], "addWorksheet")
    write_data = getExportedValue(writers[[1]], "writeData")
    save_workbook = getExportedValue(writers[[1]], "saveWorkbook")
    workbook = create_workbook()
    for (sheet in names(tables)) {
      add_worksheet(workbook, sheet)
      write_data(workbook, sheet, tables[[sheet]])
    }
    save_workbook(workbook, path, overwrite = TRUE)
    return(invisible(path))
  }
  write_xlsx = getExportedValue(writers[[2]], "write_xlsx")
  write_xlsx(tables, path)
  invisible(path)
}

.categorate_write_csv_bundle = function(tables, file_names, path, overwrite) {
  if (dir.exists(path) || file.exists(path)) {
    if (!isTRUE(overwrite)) {
      stop("Output path already exists. Use `overwrite = TRUE` to replace it.",
           call. = FALSE)
    }
    unlink(path, recursive = TRUE, force = TRUE)
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) {
    stop("Could not create output directory: ", path, call. = FALSE)
  }
  for (i in seq_along(tables)) {
    .categorate_write_csv_file(tables[[i]], file.path(path, file_names[[i]]))
  }
  invisible(path)
}

.categorate_write_csv_file = function(table, file, append = FALSE,
                                      cols = NULL,
                                      max_cell_chars = 30000) {
  if (is.null(cols)) cols = names(table)
  for (col in setdiff(cols, names(table))) table[[col]] = NA
  table = table[, cols, drop = FALSE]
  table = .categorate_export_prepare_table(table, max_cell_chars)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(table,
                     file = file,
                     sep = ",",
                     row.names = FALSE,
                     col.names = !isTRUE(append),
                     append = isTRUE(append),
                     quote = TRUE,
                     na = "",
                     qmethod = "double",
                     fileEncoding = "UTF-8")
  invisible(file)
}

.categorate_export_sheet_names = function(table_names) {
  base = vapply(table_names, function(name) {
    name = gsub("[\\[\\]\\*\\?/\\\\:]", "_", name)
    name = gsub("[[:space:]]+", "_", name)
    name = .uaf_first_non_empty_text(name, "Sheet")
    substr(name, 1, 31)
  }, character(1))

  out = character(length(base))
  seen = character()
  for (i in seq_along(base)) {
    candidate = base[[i]]
    if (candidate %in% seen) {
      suffix = 2L
      repeat {
        suffix_text = paste0("_", suffix)
        candidate = paste0(substr(base[[i]], 1, 31 - nchar(suffix_text)),
                           suffix_text)
        if (!candidate %in% seen) break
        suffix = suffix + 1L
      }
    }
    out[[i]] = candidate
    seen = c(seen, candidate)
  }
  out
}

.categorate_export_file_names = function(table_names) {
  safe = vapply(table_names, function(name) {
    name = gsub("[^A-Za-z0-9._-]+", "_", name)
    name = gsub("^_+|_+$", "", name)
    .uaf_first_non_empty_text(name, "table")
  }, character(1))
  paste0(sprintf("%02d", seq_along(safe)), "_", safe, ".csv")
}
