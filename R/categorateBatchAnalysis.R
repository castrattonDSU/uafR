#' Read a directory of resumable categorate batch files
#'
#' @description
#' `readCategorateBatchDirectory()` loads `.rds` files produced by a resumable
#' `categorate()` or categorate-enrichment batch run. It keeps successful
#' result objects and saved error objects together, then adds a batch summary so
#' downstream analysis can detect failed or incomplete batches before combining
#' tables.
#'
#' @param path Directory containing batch `.rds` files.
#' @param pattern Regular expression used to identify batch files.
#' @param min_property_ratio Minimum acceptable ratio of `PubChemProperties`
#' rows to resolved PubChem CIDs for a successful batch.
#'
#' @return A list with class `"uaf_categorate_batch_directory"` containing
#' `Batches`, `BatchSummary`, and `Source`.
#'
#' @examples
#' \dontrun{
#' batches = readCategorateBatchDirectory("categorate_batches")
#' batches$BatchSummary
#' }
#'
#' @export
readCategorateBatchDirectory = function(path,
                                        pattern = "categorate_batch_[0-9]+[.]rds$",
                                        min_property_ratio = 0.9) {
  path = .uaf_non_empty(path)
  if (length(path) != 1) {
    stop("`path` must be one batch directory.", call. = FALSE)
  }
  if (!dir.exists(path)) {
    stop("Batch directory does not exist: ", path, call. = FALSE)
  }
  files = list.files(path, pattern = pattern, full.names = TRUE)
  files = sort(files)
  if (length(files) < 1) {
    stop("No categorate batch `.rds` files were found in: ", path,
         call. = FALSE)
  }
  batches = lapply(files, readRDS)
  indices = .categorate_batch_file_indices(files)
  names(batches) = sprintf("batch_%04d", indices)
  attr(batches, "batch_files") = normalizePath(files, winslash = "/",
                                               mustWork = FALSE)
  attr(batches, "batch_indices") = indices

  out = list(
    Batches = batches,
    BatchSummary = summarizeCategorateBatches(
      batches,
      min_property_ratio = min_property_ratio
    ),
    Source = data.frame(
      Path = normalizePath(path, winslash = "/", mustWork = FALSE),
      Pattern = pattern,
      BatchCount = length(batches),
      ReadAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      stringsAsFactors = FALSE
    )
  )
  class(out) = c("uaf_categorate_batch_directory", "list")
  out
}

#' Summarize resumable categorate batches
#'
#' @description
#' `summarizeCategorateBatches()` reports which batch files are usable, which
#' contain saved errors, and whether PubChem property coverage meets the
#' requested threshold. Use this before combining large enrichment outputs.
#'
#' @param categorate_batches A batch directory path, object returned by
#' `readCategorateBatchDirectory()`, list of categorate result objects, or
#' character vector of `.rds` files.
#' @param min_property_ratio Minimum acceptable ratio of `PubChemProperties`
#' rows to resolved PubChem CIDs.
#'
#' @return A data frame with one row per batch.
#'
#' @examples
#' \dontrun{
#' summarizeCategorateBatches("categorate_batches")
#' }
#'
#' @export
summarizeCategorateBatches = function(categorate_batches,
                                      min_property_ratio = 0.9) {
  input = .categorate_batch_input(categorate_batches)
  .categorate_batch_summary_from_input(input,
                                       min_property_ratio = min_property_ratio)
}

.categorate_batch_summary_from_input = function(input,
                                                min_property_ratio = 0.9) {
  batches = input$batches
  rows = vector("list", length(batches))
  for (i in seq_along(batches)) {
    rows[[i]] = .categorate_batch_summary_row(
      batches[[i]],
      batch_index = input$indices[[i]],
      batch_name = names(batches)[[i]],
      batch_file = input$files[[i]],
      min_property_ratio = min_property_ratio
    )
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

#' Combine data frames from categorate batch results
#'
#' @description
#' `combineCategorateTables()` binds selected result tables across successful
#' categorate batches. Columns are aligned before binding, and batch metadata is
#' added by default so every row remains traceable to the source `.rds` file.
#'
#' @param categorate_batches A batch directory path, object returned by
#' `readCategorateBatchDirectory()`, list of categorate result objects, or
#' character vector of `.rds` files.
#' @param tables Character vector of categorate table names to combine. If
#' `NULL`, a curated analysis-ready set is used.
#' @param include_batch_metadata Logical. If `TRUE`, add
#' `CategorateBatchIndex`, `CategorateBatchName`, and `CategorateBatchFile`.
#' @param include_empty Logical. If `TRUE`, include empty tables when present.
#' @param min_property_ratio Minimum acceptable ratio of `PubChemProperties`
#' rows to resolved PubChem CIDs.
#' @param max_cell_chars Maximum characters retained in a single character cell.
#'
#' @return A named list of combined data frames.
#'
#' @examples
#' \dontrun{
#' combined = combineCategorateTables(
#'   "categorate_batches",
#'   tables = c("ChemicalTraits", "PubChemProperties")
#' )
#' combined$ChemicalTraits
#' }
#'
#' @export
combineCategorateTables = function(categorate_batches,
                                   tables = NULL,
                                   include_batch_metadata = TRUE,
                                   include_empty = FALSE,
                                   min_property_ratio = 0.9,
                                   max_cell_chars = 30000) {
  input = .categorate_batch_input(categorate_batches)
  batches = input$batches
  summary = .categorate_batch_summary_from_input(input,
                                                 min_property_ratio =
                                                   min_property_ratio)
  ok = summary$Status %in% c("ok", "incomplete")
  batches = batches[ok]
  if (length(batches) < 1) {
    stop("No successful categorate batches were available to combine.",
         call. = FALSE)
  }
  input$files = input$files[ok]
  input$indices = input$indices[ok]
  requested = .uaf_non_empty(tables)
  if (length(requested) < 1) {
    requested = .categorate_batch_default_tables()
  }

  out = list()
  for (table_name in requested) {
    pieces = list()
    for (i in seq_along(batches)) {
      table = batches[[i]][[table_name]]
      if (!is.data.frame(table)) next
      if (!isTRUE(include_empty) && nrow(table) < 1) next
      table = .categorate_export_prepare_table(table, max_cell_chars)
      if (isTRUE(include_batch_metadata)) {
        table = cbind(
          data.frame(
            CategorateBatchIndex = input$indices[[i]],
            CategorateBatchName = names(batches)[[i]],
            CategorateBatchFile = input$files[[i]],
            stringsAsFactors = FALSE
          ),
          table,
          stringsAsFactors = FALSE
        )
      }
      pieces[[length(pieces) + 1L]] = table
    }
    if (length(pieces) < 1) {
      if (isTRUE(include_empty)) out[[table_name]] = data.frame()
      next
    }
    out[[table_name]] = .categorate_batch_bind_rows(pieces)
  }
  out
}

#' Export a plant chemistry analysis bundle from batch outputs
#'
#' @description
#' `exportPlantChemistryAnalysisBundle()` writes a researcher-facing bundle that
#' combines resumable categorate batches with optional plant-chemistry and
#' Tanimoto tables. It is designed for large species-first projects where
#' enrichment batches, species-compound memberships, and pairwise chemistry
#' summaries need to be handed to a downstream analysis workspace.
#'
#' @param categorate_batches A batch directory path, object returned by
#' `readCategorateBatchDirectory()`, list of categorate result objects, or
#' character vector of `.rds` files.
#' @param path Output directory for `format = "csv"` or `.xlsx` path for
#' `format = "xlsx"`.
#' @param plant_membership Optional data frame or CSV path containing
#' species-compound membership rows.
#' @param species_pair_tanimoto Optional data frame or CSV path containing
#' species-pair Tanimoto summaries.
#' @param resolved_compounds Optional data frame or CSV path containing unique
#' resolved compounds.
#' @param pubchem_fingerprints Optional data frame or CSV path containing
#' PubChem fingerprint records.
#' @param file_references Optional character vector or data frame of large
#' external files to record in the bundle manifest, such as compressed
#' compound-pair Tanimoto tables.
#' @param plant_list Optional character vector, data frame, or CSV path
#' containing the full plant list used by the project. When supplied, the
#' finalized bundle includes a missing-chemistry coverage table.
#' @param metadata Optional plant metadata data frame or CSV path. If it
#' contains `species`, `accepted_species_name`, `genus`, or `family`, those
#' fields are used only for transparent taxonomy/status joins and are never
#' fabricated.
#' @param project_id Optional project label written to bundle documentation and
#' provenance text.
#' @param tables Categorate result table names to combine and export. If `NULL`,
#' a curated analysis-ready set is used.
#' @param format Export format: `"csv"`, `"xlsx"`, or `"auto"`.
#' @param include_empty Logical. If `TRUE`, export empty combined tables.
#' @param min_property_ratio Minimum acceptable ratio of `PubChemProperties`
#' rows to resolved PubChem CIDs.
#' @param overwrite Logical. If `TRUE`, replace an existing output.
#' @param max_cell_chars Maximum characters retained in a single character cell.
#' @param finalize Logical. If `TRUE` and `format = "csv"`, add enriched
#' analysis-ready tables, bundle documentation, and validation summaries.
#' @param validate_export Logical. If `TRUE`, run bundle CSV validation after
#' finalization.
#'
#' @return A manifest data frame describing exported tables.
#'
#' @examples
#' \dontrun{
#' manifest = exportPlantChemistryAnalysisBundle(
#'   categorate_batches = "categorate_batches",
#'   path = "plant_chemistry_analysis_bundle",
#'   plant_membership = "dsi_species_compound_membership.csv",
#'   species_pair_tanimoto = "dsi_species_pair_tanimoto_summary.csv",
#'   overwrite = TRUE
#' )
#' manifest
#' }
#'
#' @export
exportPlantChemistryAnalysisBundle = function(categorate_batches,
                                              path,
                                              plant_membership = NULL,
                                              species_pair_tanimoto = NULL,
                                              resolved_compounds = NULL,
                                              pubchem_fingerprints = NULL,
                                              file_references = NULL,
                                              plant_list = NULL,
                                              metadata = NULL,
                                              project_id = NULL,
                                              tables = NULL,
                                              format = c("auto", "xlsx", "csv"),
                                              include_empty = FALSE,
                                              min_property_ratio = 0.9,
                                              overwrite = FALSE,
                                              max_cell_chars = 30000,
                                              finalize = TRUE,
                                              validate_export = TRUE) {
  format = match.arg(format)
  if (missing(path) || length(.uaf_non_empty(path)) != 1) {
    stop("`path` is required.", call. = FALSE)
  }
  resolved = .categorate_export_resolve_path(path, format)
  if (identical(resolved$format, "csv") &&
      .categorate_batch_can_stream(categorate_batches)) {
    return(.categorate_analysis_export_csv_stream(
      categorate_batches = categorate_batches,
      path = resolved$path,
      plant_membership = plant_membership,
      species_pair_tanimoto = species_pair_tanimoto,
      resolved_compounds = resolved_compounds,
      pubchem_fingerprints = pubchem_fingerprints,
      file_references = file_references,
      plant_list = plant_list,
      metadata = metadata,
      project_id = project_id,
      tables = tables,
      include_empty = include_empty,
      min_property_ratio = min_property_ratio,
      overwrite = overwrite,
      max_cell_chars = max_cell_chars,
      finalize = finalize,
      validate_export = validate_export
    ))
  }
  batch_summary = summarizeCategorateBatches(
    categorate_batches,
    min_property_ratio = min_property_ratio
  )
  combined = combineCategorateTables(
    categorate_batches = categorate_batches,
    tables = tables,
    include_batch_metadata = TRUE,
    include_empty = include_empty,
    min_property_ratio = min_property_ratio,
    max_cell_chars = max_cell_chars
  )

  export_tables = list(BatchSummary = batch_summary)
  optional_tables = list(
    PlantCompoundMembership =
      .categorate_analysis_optional_table(plant_membership,
                                          "plant_membership"),
    PlantPairTanimotoSummary =
      .categorate_analysis_optional_table(species_pair_tanimoto,
                                          "species_pair_tanimoto"),
    ResolvedCompounds =
      .categorate_analysis_optional_table(resolved_compounds,
                                          "resolved_compounds"),
    PubChemFingerprints =
      .categorate_analysis_optional_table(pubchem_fingerprints,
                                          "pubchem_fingerprints"),
    FileReferences =
      .categorate_analysis_file_references(file_references)
  )
  optional_tables = optional_tables[
    vapply(optional_tables, is.data.frame, logical(1))
  ]
  optional_tables = optional_tables[
    isTRUE(include_empty) |
      vapply(optional_tables, nrow, integer(1)) > 0
  ]
  export_tables = c(export_tables, optional_tables, combined)

  if (!isTRUE(include_empty)) {
    export_tables = export_tables[
      vapply(export_tables, nrow, integer(1)) > 0
    ]
  }
  if (length(export_tables) < 1) {
    stop("No exportable analysis tables were found.", call. = FALSE)
  }

  all_table_names = c("ExportManifest", names(export_tables))
  sheet_names = .categorate_export_sheet_names(all_table_names)
  file_names = .categorate_export_file_names(all_table_names)
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  manifest = data.frame(
    Table = all_table_names,
    SheetName = sheet_names,
    FileName = if (resolved$format == "csv") file_names else NA_character_,
    RowCount = c(NA_integer_, vapply(export_tables, nrow, integer(1))),
    ColumnCount = c(NA_integer_, vapply(export_tables, ncol, integer(1))),
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
    if (isTRUE(finalize)) {
      manifest = finalizePlantChemistryAnalysisBundle(
        path = resolved$path,
        plant_list = plant_list,
        metadata = metadata,
        project_id = project_id,
        overwrite = TRUE,
        validate_export = validate_export,
        max_cell_chars = max_cell_chars
      )
    }
  }
  row.names(manifest) = NULL
  manifest
}

.categorate_batch_input = function(x) {
  if (inherits(x, "uaf_categorate_batch_directory")) {
    batches = x$Batches
    files = attr(batches, "batch_files")
    indices = attr(batches, "batch_indices")
    return(.categorate_batch_standardize_input(batches, files, indices))
  }
  if (is.list(x) && !is.data.frame(x) && "Batches" %in% names(x) &&
      is.list(x$Batches)) {
    batches = x$Batches
    files = attr(batches, "batch_files")
    indices = attr(batches, "batch_indices")
    return(.categorate_batch_standardize_input(batches, files, indices))
  }
  if (is.character(x)) {
    x = .uaf_non_empty(x)
    if (length(x) == 1 && dir.exists(x)) {
      files = sort(list.files(x, pattern = "[.]rds$", full.names = TRUE))
    } else {
      files = x
    }
    if (length(files) < 1 || any(!file.exists(files))) {
      stop("All batch file paths must exist.", call. = FALSE)
    }
    files = sort(files)
    batches = lapply(files, readRDS)
    indices = .categorate_batch_file_indices(files)
    names(batches) = sprintf("batch_%04d", indices)
    files = normalizePath(files, winslash = "/", mustWork = FALSE)
    return(.categorate_batch_standardize_input(batches, files, indices))
  }
  if (is.list(x) && !is.data.frame(x)) {
    if (.categorate_batch_is_result(x)) {
      return(.categorate_batch_standardize_input(list(batch_0001 = x),
                                                 NA_character_, 1L))
    }
    if (all(vapply(x, .categorate_batch_is_result, logical(1)))) {
      return(.categorate_batch_standardize_input(
        x,
        attr(x, "batch_files"),
        attr(x, "batch_indices")
      ))
    }
  }
  stop("`categorate_batches` must be a batch directory, `.rds` file vector, ",
       "categorate result, or list of categorate results.", call. = FALSE)
}

.categorate_batch_standardize_input = function(batches, files, indices) {
  n = length(batches)
  if (n < 1) stop("At least one categorate batch is required.", call. = FALSE)
  if (is.null(indices) || length(indices) != n) indices = seq_len(n)
  indices = suppressWarnings(as.integer(indices))
  indices[is.na(indices)] = seq_len(n)[is.na(indices)]
  if (is.null(files) || length(files) != n) files = rep(NA_character_, n)
  if (is.null(names(batches)) || any(names(batches) == "")) {
    names(batches) = sprintf("batch_%04d", indices)
  }
  list(batches = batches, files = files, indices = indices)
}

.categorate_batch_is_result = function(x) {
  inherits(x, "error") ||
    (is.list(x) && !is.data.frame(x) &&
       any(c("PubChemIdentity", "ChemicalTraits", "ValidationSummary",
             "ChemicalTraitMatrix", "SourceCoverage") %in% names(x)))
}

.categorate_batch_file_indices = function(files) {
  base = basename(files)
  hits = regmatches(base, regexpr("[0-9]+(?=[.]rds$)", base, perl = TRUE))
  out = suppressWarnings(as.integer(hits))
  bad = is.na(out)
  out[bad] = seq_along(files)[bad]
  out
}

.categorate_batch_summary_row = function(result, batch_index, batch_name,
                                         batch_file, min_property_ratio) {
  if (inherits(result, "error")) {
    return(data.frame(
      BatchIndex = batch_index,
      BatchName = batch_name,
      BatchFile = batch_file,
      Status = "error",
      ValidationStatus = NA_character_,
      QueryCount = NA_integer_,
      ResolvedCIDCount = NA_integer_,
      PubChemPropertyRows = NA_integer_,
      PubChemPropertyRatio = NA_real_,
      ChemicalTraitRows = NA_integer_,
      ChemicalTraitMatrixRows = NA_integer_,
      ValidationIssueCount = NA_integer_,
      QualityIssue = conditionMessage(result),
      Error = conditionMessage(result),
      stringsAsFactors = FALSE
    ))
  }

  identity = result$PubChemIdentity
  properties = result$PubChemProperties
  resolved = if (is.data.frame(identity) && "CID" %in% names(identity)) {
    sum(!is.na(identity$CID))
  } else {
    0L
  }
  property_rows = if (is.data.frame(properties)) nrow(properties) else 0L
  property_ratio = if (resolved > 0) property_rows / resolved else 1
  incomplete = resolved > 0 && property_ratio < min_property_ratio
  validation = .categorate_batch_validation_status(result)
  issue_count = if (is.data.frame(result$ValidationIssues)) {
    nrow(result$ValidationIssues)
  } else {
    NA_integer_
  }
  query_count = if (is.data.frame(identity) && "Query" %in% names(identity)) {
    length(unique(.uaf_non_empty(identity$Query)))
  } else {
    NA_integer_
  }

  data.frame(
    BatchIndex = batch_index,
    BatchName = batch_name,
    BatchFile = batch_file,
    Status = if (incomplete) "incomplete" else "ok",
    ValidationStatus = validation,
    QueryCount = query_count,
    ResolvedCIDCount = resolved,
    PubChemPropertyRows = property_rows,
    PubChemPropertyRatio = round(property_ratio, 4),
    ChemicalTraitRows = .categorate_batch_table_rows(result,
                                                     "ChemicalTraits"),
    ChemicalTraitMatrixRows = .categorate_batch_table_rows(
      result,
      "ChemicalTraitMatrix"
    ),
    ValidationIssueCount = issue_count,
    QualityIssue = if (incomplete) {
      paste0("Only ", property_rows, " PubChem property rows for ",
             resolved, " resolved CIDs.")
    } else {
      NA_character_
    },
    Error = NA_character_,
    stringsAsFactors = FALSE
  )
}

.categorate_batch_validation_status = function(result) {
  if (is.data.frame(result$ValidationSummary) &&
      "Status" %in% names(result$ValidationSummary) &&
      nrow(result$ValidationSummary) > 0) {
    return(paste(unique(.uaf_non_empty(result$ValidationSummary$Status)),
                 collapse = "; "))
  }
  NA_character_
}

.categorate_batch_table_rows = function(result, table_name) {
  table = result[[table_name]]
  if (is.data.frame(table)) return(nrow(table))
  NA_integer_
}

.categorate_batch_default_tables = function() {
  c("PubChemIdentity", "PubChemProperties", "ChemicalTraitReport",
    "DerivedGroups", "ChemicalTraitMatrix", "ChemicalTraitOntologyMatrix",
    "ChemicalTraitSummary", "ChemicalTraitSimilarity", "ChemicalTraits",
    "ChemicalTraitOntology", "ChemicalTraitEvidence", "ChemicalMeasurements",
    "ChemicalMeasurementSummary", "ChemicalHazards", "ChemicalUses",
    "ChemicalBioassays", "ChemicalBioactivities", "ChemicalTargets",
    "ChemicalPotencies", "ChemicalTaxonomy", "ChemicalOccurrences",
    "ChemicalPathwayRoles", "KEGGReactionParticipants", "ChemicalClasses",
    "SourceCoverage", "TableQuality", "SourceDiagnostics",
    "ValidationSummary", "ValidationIssues", "DataDictionary", "Provenance")
}

.categorate_batch_bind_rows = function(rows) {
  cols = unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned = lapply(rows, function(row) {
    for (col in setdiff(cols, names(row))) row[[col]] = NA
    row[, cols, drop = FALSE]
  })
  out = do.call(rbind, aligned)
  row.names(out) = NULL
  out
}

.categorate_analysis_optional_table = function(x, label) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  path = .uaf_non_empty(x)
  if (length(path) != 1) {
    stop("`", label, "` must be a data frame or one file path.",
         call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input file for `", label, "` does not exist: ", path,
         call. = FALSE)
  }
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    out = readRDS(path)
    if (!is.data.frame(out)) {
      stop("RDS input for `", label, "` must contain a data frame.",
           call. = FALSE)
    }
    return(out)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.categorate_analysis_file_references = function(x) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(x)
  original_names = names(x)
  paths = .uaf_non_empty(x)
  if (length(paths) < 1) return(NULL)
  names_in = original_names[match(paths, as.character(x))]
  data.frame(
    Reference = if (!is.null(names_in) && length(names_in) == length(paths)) {
      ifelse(names_in == "", basename(paths), names_in)
    } else {
      basename(paths)
    },
    Path = normalizePath(paths, winslash = "/", mustWork = FALSE),
    Exists = file.exists(paths),
    SizeBytes = ifelse(file.exists(paths), file.info(paths)$size, NA_real_),
    stringsAsFactors = FALSE
  )
}

.categorate_batch_can_stream = function(x) {
  !is.null(.categorate_batch_file_input(x, error = FALSE))
}

.categorate_batch_file_input = function(x, error = TRUE) {
  if (!is.character(x)) {
    if (isTRUE(error)) {
      stop("Streaming batch export requires a batch directory or `.rds` file ",
           "paths.", call. = FALSE)
    }
    return(NULL)
  }
  x = .uaf_non_empty(x)
  if (length(x) == 1 && dir.exists(x)) {
    files = sort(list.files(x, pattern = "[.]rds$", full.names = TRUE))
  } else {
    files = sort(x)
  }
  if (length(files) < 1 || any(!file.exists(files))) {
    if (isTRUE(error)) stop("All batch file paths must exist.", call. = FALSE)
    return(NULL)
  }
  files = normalizePath(files, winslash = "/", mustWork = FALSE)
  indices = .categorate_batch_file_indices(files)
  list(files = files,
       indices = indices,
       names = sprintf("batch_%04d", indices))
}

.categorate_analysis_export_csv_stream = function(categorate_batches,
                                                  path,
                                                  plant_membership,
                                                  species_pair_tanimoto,
                                                  resolved_compounds,
                                                  pubchem_fingerprints,
                                                  file_references,
                                                  plant_list,
                                                  metadata,
                                                  project_id,
                                                  tables,
                                                  include_empty,
                                                  min_property_ratio,
                                                  overwrite,
                                                  max_cell_chars,
                                                  finalize,
                                                  validate_export) {
  file_input = .categorate_batch_file_input(categorate_batches)
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

  batch_summary = .categorate_batch_summary_from_files(
    file_input,
    min_property_ratio = min_property_ratio
  )
  optional_tables = list(
    PlantCompoundMembership =
      .categorate_analysis_optional_table(plant_membership,
                                          "plant_membership"),
    PlantPairTanimotoSummary =
      .categorate_analysis_optional_table(species_pair_tanimoto,
                                          "species_pair_tanimoto"),
    ResolvedCompounds =
      .categorate_analysis_optional_table(resolved_compounds,
                                          "resolved_compounds"),
    PubChemFingerprints =
      .categorate_analysis_optional_table(pubchem_fingerprints,
                                          "pubchem_fingerprints"),
    FileReferences =
      .categorate_analysis_file_references(file_references)
  )
  optional_tables = optional_tables[
    vapply(optional_tables, is.data.frame, logical(1))
  ]
  optional_tables = optional_tables[
    isTRUE(include_empty) |
      vapply(optional_tables, nrow, integer(1)) > 0
  ]

  requested = .uaf_non_empty(tables)
  if (length(requested) < 1) requested = .categorate_batch_default_tables()
  stream_plan = .categorate_batch_stream_plan(
    file_input = file_input,
    batch_summary = batch_summary,
    tables = requested,
    include_empty = include_empty,
    max_cell_chars = max_cell_chars
  )
  stream_tables = names(stream_plan$row_counts)[
    stream_plan$row_counts > 0 | isTRUE(include_empty)
  ]

  fixed_tables = c(list(BatchSummary = batch_summary), optional_tables)
  fixed_tables = fixed_tables[
    isTRUE(include_empty) | vapply(fixed_tables, nrow, integer(1)) > 0
  ]
  table_names = c(names(fixed_tables), stream_tables)
  all_table_names = c("ExportManifest", table_names)
  sheet_names = .categorate_export_sheet_names(all_table_names)
  file_names = .categorate_export_file_names(all_table_names)
  names(file_names) = all_table_names

  for (table_name in names(fixed_tables)) {
    table = .categorate_export_prepare_table(fixed_tables[[table_name]],
                                             max_cell_chars)
    .categorate_stream_write_csv(
      table,
      file.path(path, file_names[[table_name]]),
      append = FALSE,
      cols = names(table)
    )
  }
  .categorate_batch_stream_write_tables(
    file_input = file_input,
    batch_summary = batch_summary,
    tables = stream_tables,
    plan = stream_plan,
    path = path,
    file_names = file_names,
    max_cell_chars = max_cell_chars
  )

  row_counts = c(vapply(fixed_tables, nrow, integer(1)),
                 stream_plan$row_counts[stream_tables])
  col_counts = c(vapply(fixed_tables, ncol, integer(1)),
                 vapply(stream_plan$cols[stream_tables], length, integer(1)))
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  manifest = data.frame(
    Table = all_table_names,
    SheetName = sheet_names,
    FileName = file_names,
    RowCount = c(NA_integer_, as.integer(row_counts)),
    ColumnCount = c(NA_integer_, as.integer(col_counts)),
    Format = "csv",
    OutputPath = path,
    CreatedAt = created_at,
    stringsAsFactors = FALSE
  )
  manifest$RowCount[manifest$Table == "ExportManifest"] = nrow(manifest)
  manifest$ColumnCount[manifest$Table == "ExportManifest"] = ncol(manifest)
  .categorate_stream_write_csv(
    manifest,
    file.path(path, file_names[["ExportManifest"]]),
    append = FALSE,
    cols = names(manifest)
  )
  if (isTRUE(finalize)) {
    manifest = finalizePlantChemistryAnalysisBundle(
      path = path,
      plant_list = plant_list,
      metadata = metadata,
      project_id = project_id,
      overwrite = TRUE,
      validate_export = validate_export,
      max_cell_chars = max_cell_chars
    )
  }
  row.names(manifest) = NULL
  manifest
}

.categorate_batch_summary_from_files = function(file_input,
                                                min_property_ratio = 0.9) {
  rows = vector("list", length(file_input$files))
  for (i in seq_along(file_input$files)) {
    result = readRDS(file_input$files[[i]])
    rows[[i]] = .categorate_batch_summary_row(
      result,
      batch_index = file_input$indices[[i]],
      batch_name = file_input$names[[i]],
      batch_file = file_input$files[[i]],
      min_property_ratio = min_property_ratio
    )
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.categorate_batch_stream_plan = function(file_input, batch_summary, tables,
                                         include_empty, max_cell_chars) {
  cols = stats::setNames(vector("list", length(tables)), tables)
  row_counts = stats::setNames(rep(0L, length(tables)), tables)
  ok = batch_summary$Status %in% c("ok", "incomplete")
  for (i in which(ok)) {
    result = readRDS(file_input$files[[i]])
    if (inherits(result, "error")) next
    for (table_name in tables) {
      table = result[[table_name]]
      if (!is.data.frame(table)) next
      if (!isTRUE(include_empty) && nrow(table) < 1) next
      table = .categorate_export_prepare_table(table, max_cell_chars)
      table = .categorate_batch_add_metadata(
        table,
        batch_index = file_input$indices[[i]],
        batch_name = file_input$names[[i]],
        batch_file = file_input$files[[i]]
      )
      cols[[table_name]] = unique(c(cols[[table_name]], names(table)))
      row_counts[[table_name]] = row_counts[[table_name]] + nrow(table)
    }
  }
  list(cols = cols, row_counts = row_counts)
}

.categorate_batch_stream_write_tables = function(file_input, batch_summary,
                                                 tables, plan, path,
                                                 file_names,
                                                 max_cell_chars) {
  if (length(tables) < 1) return(invisible(FALSE))
  ok = batch_summary$Status %in% c("ok", "incomplete")
  for (table_name in tables) {
    out_file = file.path(path, file_names[[table_name]])
    wrote = FALSE
    for (i in which(ok)) {
      result = readRDS(file_input$files[[i]])
      if (inherits(result, "error")) next
      table = result[[table_name]]
      if (!is.data.frame(table) || nrow(table) < 1) next
      table = .categorate_export_prepare_table(table, max_cell_chars)
      table = .categorate_batch_add_metadata(
        table,
        batch_index = file_input$indices[[i]],
        batch_name = file_input$names[[i]],
        batch_file = file_input$files[[i]]
      )
      .categorate_stream_write_csv(table, out_file,
                                   append = wrote,
                                   cols = plan$cols[[table_name]])
      wrote = TRUE
    }
    if (!wrote && length(plan$cols[[table_name]]) > 0) {
      empty = as.data.frame(stats::setNames(
        rep(list(character()), length(plan$cols[[table_name]])),
        plan$cols[[table_name]]
      ))
      .categorate_stream_write_csv(empty, out_file,
                                   append = FALSE,
                                   cols = plan$cols[[table_name]])
    }
  }
  invisible(TRUE)
}

.categorate_batch_add_metadata = function(table, batch_index, batch_name,
                                          batch_file) {
  cbind(
    data.frame(
      CategorateBatchIndex = batch_index,
      CategorateBatchName = batch_name,
      CategorateBatchFile = batch_file,
      stringsAsFactors = FALSE
    ),
    table,
    stringsAsFactors = FALSE
  )
}

.categorate_stream_write_csv = function(table, file, append, cols) {
  .categorate_write_csv_file(table, file, append = append, cols = cols)
  invisible(file)
}
