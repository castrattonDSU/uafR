#' Finalize a plant chemistry analysis bundle
#'
#' @description
#' `finalizePlantChemistryAnalysisBundle()` adds the analysis-ready layer to a
#' CSV bundle written by `exportPlantChemistryAnalysisBundle()`. It does not
#' query public services. Instead, it uses the exported plant-compound
#' membership, Tanimoto summaries, PubChem fingerprints/properties, and
#' categorate-derived tables already present in the bundle.
#'
#' The finalizer writes enriched plant-compound membership, species chemistry
#' summaries, missing-species coverage, source coverage summaries, validation
#' overviews, a data dictionary, and README/methods text. It also refreshes the
#' export manifest so downstream projects can audit row and column counts.
#'
#' @param path CSV bundle directory.
#' @param plant_list Optional character vector, data frame, or CSV path
#' containing the complete project plant list. A `species` column is preferred
#' when a data frame is supplied.
#' @param metadata Optional plant metadata data frame or CSV path with a
#' `species` column. `accepted_species_name`, `genus`, and `family` are joined
#' when available.
#' @param project_id Optional project label written to documentation.
#' @param overwrite Logical. If `FALSE`, existing finalization outputs are not
#' overwritten.
#' @param validate_export Logical. If `TRUE`, run bundle validation and include
#' validation status in `12b_ValidationOverview.csv`.
#' @param max_cell_chars Maximum characters retained in a single exported cell.
#'
#' @return Updated export manifest.
#'
#' @export
finalizePlantChemistryAnalysisBundle = function(path,
                                                plant_list = NULL,
                                                metadata = NULL,
                                                project_id = NULL,
                                                overwrite = TRUE,
                                                validate_export = TRUE,
                                                max_cell_chars = 30000) {
  path = .uaf_non_empty(path)
  if (length(path) != 1 || !dir.exists(path)) {
    stop("`path` must be one existing CSV bundle directory.", call. = FALSE)
  }
  path = normalizePath(path, winslash = "/", mustWork = FALSE)
  manifest = .bundle_read_manifest(path)
  metadata = .bundle_optional_frame(metadata, "metadata")
  plant_list = .bundle_standardize_plant_list(plant_list, metadata)

  membership = .bundle_read_table(path, manifest, "PlantCompoundMembership",
                                  required = TRUE)
  resolved = .bundle_read_table(path, manifest, "ResolvedCompounds")
  fingerprints = .bundle_read_table(path, manifest, "PubChemFingerprints")
  pair_tanimoto = .bundle_read_table(path, manifest,
                                     "PlantPairTanimotoSummary")
  derived = .bundle_read_table(path, manifest, "DerivedGroups")
  properties = .bundle_read_table(path, manifest, "PubChemProperties")
  coverage = .bundle_read_table(path, manifest, "SourceCoverage")
  validation_summary = .bundle_read_table(path, manifest,
                                          "ValidationSummary")
  validation_issues = .bundle_read_table(path, manifest, "ValidationIssues")

  membership = .bundle_apply_taxonomy_status(membership, metadata)
  .bundle_write_table(path, "03_PlantCompoundMembership.csv", membership,
                      overwrite = overwrite,
                      max_cell_chars = max_cell_chars)

  enriched = .bundle_enrich_membership(
    membership = membership,
    resolved = resolved,
    fingerprints = fingerprints,
    properties = properties,
    derived = derived,
    metadata = metadata
  )
  .bundle_write_table(path, "03b_PlantCompoundMembershipEnriched.csv",
                      enriched, overwrite = overwrite,
                      max_cell_chars = max_cell_chars)

  if (is.data.frame(pair_tanimoto) && nrow(pair_tanimoto) > 0) {
    pair_tanimoto = .bundle_enhance_pair_tanimoto(pair_tanimoto)
    .bundle_write_table(path, "04_PlantPairTanimotoSummary.csv",
                        pair_tanimoto, overwrite = TRUE,
                        max_cell_chars = max_cell_chars)
  }

  source_summary = .bundle_source_coverage_summary(coverage)
  .bundle_write_table(path, "11b_SourceCoverageSummary.csv", source_summary,
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  species_summary = .bundle_species_summary(enriched)
  .bundle_write_table(path, "14_PlantChemistrySummary.csv", species_summary,
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  missing = .bundle_missing_species(plant_list, enriched, metadata)
  .bundle_write_table(path, "15_PlantChemistryMissingSpecies.csv", missing,
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  manifest = .bundle_refresh_manifest(path, project_id = project_id)
  dictionary = .bundle_data_dictionary(path, manifest)
  .bundle_write_table(path, "00_DataDictionary.csv", dictionary,
                      overwrite = TRUE, max_cell_chars = max_cell_chars)

  .bundle_write_readme(path, project_id = project_id)
  .bundle_write_methods(path, project_id = project_id)

  manifest = .bundle_refresh_manifest(path, project_id = project_id)
  validation = if (isTRUE(validate_export)) {
    validatePlantChemistryAnalysisBundle(path)
  } else {
    NULL
  }
  overview = .bundle_validation_overview(
    validation = validation,
    batch_validation = validation_summary,
    validation_issues = validation_issues
  )
  .bundle_write_table(path, "12b_ValidationOverview.csv", overview,
                      overwrite = TRUE, max_cell_chars = max_cell_chars)
  manifest = .bundle_refresh_manifest(path, project_id = project_id)
  attr(manifest, "Validation") = if (isTRUE(validate_export)) {
    validatePlantChemistryAnalysisBundle(path)
  } else {
    NULL
  }
  row.names(manifest) = NULL
  manifest
}

#' Validate a plant chemistry analysis bundle
#'
#' @description
#' `validatePlantChemistryAnalysisBundle()` checks a CSV bundle for parser
#' consistency, manifest row/column agreement, accidental row-index columns,
#' and required columns in the main analysis tables. Python `csv.reader` and
#' pandas checks are optional so package tests do not depend on a Python
#' installation, but the function records whether those checks were run.
#'
#' @param path CSV bundle directory.
#' @param use_python Logical. If `TRUE`, also validate every CSV with Python's
#' standard `csv.reader` when `python3` is available.
#' @param use_pandas Logical. If `TRUE`, also validate every CSV with
#' `pandas.read_csv()` when Python and pandas are available.
#'
#' @return A list with `Summary`, `CSVValidation`, `RequiredColumns`, and
#' `Manifest` tables.
#'
#' @export
validatePlantChemistryAnalysisBundle = function(path,
                                                use_python = FALSE,
                                                use_pandas = FALSE) {
  path = .uaf_non_empty(path)
  if (length(path) != 1 || !dir.exists(path)) {
    stop("`path` must be one existing CSV bundle directory.", call. = FALSE)
  }
  path = normalizePath(path, winslash = "/", mustWork = FALSE)
  manifest = .bundle_read_manifest(path, required = FALSE)
  csv_files = sort(list.files(path, pattern = "[.]csv$", full.names = FALSE))
  rows = lapply(csv_files, function(file_name) {
    .bundle_validate_csv_file(file.path(path, file_name), manifest,
                              use_python = use_python,
                              use_pandas = use_pandas)
  })
  csv_validation = .bundle_bind(rows)
  required = .bundle_validate_required_columns(path, manifest)
  blocking = any(csv_validation$Status == "fail") ||
    any(required$Status == "fail")
  warnings = any(csv_validation$Status == "warn") ||
    any(required$Status == "warn")
  summary = data.frame(
    CSVFileCount = length(csv_files),
    CSVPassCount = sum(csv_validation$Status == "pass", na.rm = TRUE),
    CSVWarnCount = sum(csv_validation$Status == "warn", na.rm = TRUE),
    CSVFailCount = sum(csv_validation$Status == "fail", na.rm = TRUE),
    RequiredColumnFailCount = sum(required$Status == "fail", na.rm = TRUE),
    PythonCsvChecked = .bundle_yes_no(isTRUE(use_python)),
    PandasChecked = .bundle_yes_no(isTRUE(use_pandas)),
    ExportReadyStatus = if (blocking) {
      "fail"
    } else if (warnings) {
      "warn"
    } else {
      "pass"
    },
    stringsAsFactors = FALSE
  )
  list(
    Summary = summary,
    CSVValidation = csv_validation,
    RequiredColumns = required,
    Manifest = manifest
  )
}

.bundle_read_manifest = function(path, required = TRUE) {
  file = file.path(path, "01_ExportManifest.csv")
  if (!file.exists(file)) {
    if (isTRUE(required)) {
      stop("Bundle manifest not found: ", file, call. = FALSE)
    }
    return(data.frame())
  }
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

.bundle_optional_frame = function(x, label) {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  path = .uaf_non_empty(x)
  if (length(path) != 1) {
    stop("`", label, "` must be a data frame or one file path.",
         call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("Input file for `", label, "` does not exist: ", path,
         call. = FALSE)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

.bundle_standardize_plant_list = function(plant_list, metadata) {
  if (is.null(plant_list)) {
    if (is.data.frame(metadata) && "species" %in% names(metadata)) {
      plant_list = metadata["species"]
    } else {
      return(data.frame(species = character(), stringsAsFactors = FALSE))
    }
  }
  if (is.character(plant_list) && length(plant_list) == 1 &&
      file.exists(plant_list)) {
    plant_list = utils::read.csv(plant_list, stringsAsFactors = FALSE,
                                 check.names = FALSE)
  }
  if (is.character(plant_list)) {
    out = data.frame(species = plant_list, stringsAsFactors = FALSE)
  } else if (is.data.frame(plant_list)) {
    out = as.data.frame(plant_list, stringsAsFactors = FALSE)
    if (!"species" %in% names(out)) {
      names(out)[[1]] = "species"
    }
  } else {
    stop("`plant_list` must be a character vector, data frame, or CSV path.",
         call. = FALSE)
  }
  out$species = .bundle_squish(out$species)
  out = out[!is.na(out$species) & out$species != "", , drop = FALSE]
  out = out[!duplicated(out$species), , drop = FALSE]
  row.names(out) = NULL
  out
}

.bundle_read_table = function(path, manifest, table, required = FALSE) {
  if (!is.data.frame(manifest) || nrow(manifest) < 1 ||
      !"Table" %in% names(manifest) || !"FileName" %in% names(manifest)) {
    if (isTRUE(required)) stop("Manifest is missing or malformed.",
                               call. = FALSE)
    return(data.frame())
  }
  hit = which(manifest$Table == table)
  if (length(hit) < 1) {
    if (isTRUE(required)) stop("Required bundle table not found: ", table,
                               call. = FALSE)
    return(data.frame())
  }
  file = file.path(path, manifest$FileName[[hit[[1]]]])
  if (!file.exists(file)) {
    if (isTRUE(required)) stop("Required bundle file not found: ", file,
                               call. = FALSE)
    return(data.frame())
  }
  utils::read.csv(file, stringsAsFactors = FALSE, check.names = FALSE)
}

.bundle_write_table = function(path, file_name, table, overwrite = TRUE,
                               max_cell_chars = 30000) {
  file = file.path(path, file_name)
  if (file.exists(file) && !isTRUE(overwrite)) {
    stop("Output file exists. Use `overwrite = TRUE`: ", file,
         call. = FALSE)
  }
  .categorate_write_csv_file(table, file, max_cell_chars = max_cell_chars)
}

.bundle_apply_taxonomy_status = function(membership, metadata = NULL) {
  membership = as.data.frame(membership, stringsAsFactors = FALSE)
  if (!"species" %in% names(membership)) {
    stop("PlantCompoundMembership must contain `species`.", call. = FALSE)
  }
  if (!"genus" %in% names(membership)) {
    membership$genus = .bundle_genus(membership$species)
  }
  if (!"family" %in% names(membership)) membership$family = NA_character_
  if (!"source_database" %in% names(membership)) {
    membership$source_database = "unknown"
  }
  membership$source_database = .bundle_first_non_empty(
    membership$source_database, "unknown"
  )
  if (!"evidence_tier" %in% names(membership)) {
    membership$evidence_tier = "unresolved"
  }
  membership$evidence_tier = .bundle_first_non_empty(
    membership$evidence_tier, "unresolved"
  )
  if (is.data.frame(metadata) && nrow(metadata) > 0 &&
      "species" %in% names(metadata)) {
    meta = metadata[!duplicated(metadata$species), , drop = FALSE]
    idx = match(membership$species, meta$species)
    if ("accepted_species_name" %in% names(meta) &&
        !"accepted_species_name" %in% names(membership)) {
      membership$accepted_species_name = .bundle_index(meta, idx,
                                                       "accepted_species_name")
    }
    if ("genus" %in% names(meta)) {
      membership$genus = .bundle_first_non_empty(
        membership$genus,
        .bundle_index(meta, idx, "genus")
      )
    }
    if ("family" %in% names(meta)) {
      membership$family = .bundle_first_non_empty(
        membership$family,
        .bundle_index(meta, idx, "family")
      )
    }
  }
  family_has_value = length(.uaf_non_empty(membership$family)) > 0
  membership$taxonomy_family_status = if (family_has_value) {
    ifelse(is.na(membership$family) | membership$family == "",
           "family_missing_for_record", "family_supplied")
  } else {
    "family_not_available_from_supplied_sources"
  }
  membership
}

.bundle_enrich_membership = function(membership, resolved, fingerprints,
                                     properties, derived, metadata = NULL) {
  out = .bundle_apply_taxonomy_status(membership, metadata)
  out$row_key = seq_len(nrow(out))
  for (col in c("compound_id", "compound_name", "compound_name_clean",
                "CID", "InChIKey", "SMILES", "MolecularFormula")) {
    if (!col %in% names(out)) out[[col]] = NA_character_
  }

  resolved = .bundle_unique_by(resolved, "compound_id")
  fingerprints = .bundle_unique_by(fingerprints, "compound_id")
  if (nrow(resolved) > 0) {
    idx = match(out$compound_id, resolved$compound_id)
    out$compound_name = .bundle_first_non_empty(
      out$compound_name, .bundle_index(resolved, idx, "compound_name")
    )
    out$compound_name_clean = .bundle_first_non_empty(
      out$compound_name_clean,
      .bundle_index(resolved, idx, "compound_name_clean")
    )
    out$InChIKey = .bundle_first_non_empty(
      out$InChIKey, .bundle_index(resolved, idx, "InChIKey")
    )
    out$SMILES = .bundle_first_non_empty(
      out$SMILES, .bundle_index(resolved, idx, "SMILES")
    )
    out$MolecularFormula = .bundle_first_non_empty(
      out$MolecularFormula, .bundle_index(resolved, idx, "MolecularFormula")
    )
    out$resolution_source = .bundle_first_non_empty(
      .bundle_col_or(out, "resolution_source", NA_character_),
      .bundle_index(resolved, idx, "resolution_source")
    )
    out$identity_review_required = .bundle_first_non_empty(
      .bundle_col_or(out, "identity_review_required", NA_character_),
      .bundle_index(resolved, idx, "identity_review_required")
    )
  }

  if (nrow(fingerprints) > 0) {
    idx = match(out$compound_id, fingerprints$compound_id)
    out$pubchem_cid = .bundle_index(fingerprints, idx, "pubchem_cid")
    out$Fingerprint2D = .bundle_index(fingerprints, idx, "Fingerprint2D")
    out$CID = .bundle_first_non_empty(out$CID, out$pubchem_cid,
                                      .bundle_index(fingerprints, idx, "CID"))
    out$InChIKey = .bundle_first_non_empty(
      out$InChIKey, .bundle_index(fingerprints, idx, "InChIKey"),
      .bundle_index(fingerprints, idx, "InChIKey_source")
    )
    out$SMILES = .bundle_first_non_empty(
      out$SMILES,
      .bundle_index(fingerprints, idx, "IsomericSMILES"),
      .bundle_index(fingerprints, idx, "CanonicalSMILES"),
      .bundle_index(fingerprints, idx, "SMILES")
    )
    out$MolecularFormula = .bundle_first_non_empty(
      out$MolecularFormula,
      .bundle_index(fingerprints, idx, "MolecularFormula"),
      .bundle_index(fingerprints, idx, "MolecularFormula_source")
    )
  } else {
    out$pubchem_cid = NA_character_
    out$Fingerprint2D = NA_character_
  }
  out$has_fingerprint = .bundle_yes_no(!is.na(out$Fingerprint2D) &
                                         out$Fingerprint2D != "")

  properties = .bundle_prepare_cid_table(properties)
  derived = .bundle_prepare_cid_table(derived)
  prop_idx = match(.bundle_cid_key(out), properties$cid_key)
  der_idx = match(.bundle_cid_key(out), derived$cid_key)
  prop_cols = c("MolecularWeight", "ExactMass", "XLogP", "TPSA",
                "HBondDonorCount", "HBondAcceptorCount",
                "RotatableBondCount", "HeavyAtomCount")
  for (col in prop_cols) out[[col]] = .bundle_index(properties, prop_idx, col)
  out$MolecularFormula = .bundle_first_non_empty(
    out$MolecularFormula, .bundle_index(properties, prop_idx,
                                        "MolecularFormula")
  )
  out$InChIKey = .bundle_first_non_empty(
    out$InChIKey, .bundle_index(properties, prop_idx, "InChIKey")
  )
  out$SMILES = .bundle_first_non_empty(
    out$SMILES, .bundle_index(properties, prop_idx, "IsomericSMILES"),
    .bundle_index(properties, prop_idx, "CanonicalSMILES"),
    .bundle_index(properties, prop_idx, "SMILES")
  )

  derived_cols = c(
    "is_natural_product", "natural_product_superclasses",
    "natural_product_classes", "natural_product_subclasses",
    "occurrence_count", "organism_count", "kingdom_count", "family_count",
    "genus_count", "kingdoms_observed", "families_observed",
    "genera_observed", "dominant_kingdom", "dominant_family",
    "taxonomic_breadth", "is_plant_occurring", "is_fungal_occurring",
    "is_bacterial_occurring", "has_kegg_pathway", "has_bioactivity",
    "has_active_bioactivity", "has_safety_hazard", "has_literature",
    "pathway_count", "source_count", "molecular_size_bin", "polarity_bin",
    "lipophilicity_bin", "volatility_proxy", "oxygenated", "nitrogenous",
    "sulfur_containing", "halogenated", "metabolic_context",
    "biomedical_context", "ecological_context", "sensory_context",
    "safety_context", "analytical_context", "bioactivity_context",
    "confidence_score", "evidence_score"
  )
  for (col in derived_cols) out[[col]] = .bundle_index(derived, der_idx, col)

  out$context_known_record = .bundle_yes_no(
    .bundle_known(.bundle_col_or(out, "plant_part_group", NA_character_)) |
      .bundle_known(.bundle_col_or(out, "tissue_group", NA_character_)) |
      .bundle_known(.bundle_col_or(out, "method_group", NA_character_))
  )
  out$plant_part_known = .bundle_yes_no(.bundle_known(
    .bundle_col_or(out, "plant_part_group", NA_character_)
  ))
  out$tissue_known = .bundle_yes_no(.bundle_known(
    .bundle_col_or(out, "tissue_group", NA_character_)
  ))
  out$method_known = .bundle_yes_no(.bundle_known(
    .bundle_col_or(out, "method_group", NA_character_)
  ))

  comparability = .bundle_membership_comparability(out, derived)
  if (nrow(comparability) > 0) {
    key = .bundle_occurrence_key(out)
    comp_key = .bundle_occurrence_key(comparability)
    idx = match(key, comp_key)
    for (col in c("metabolism_domain", "biosynthetic_family",
                  "chemical_behavior", "comparison_scope",
                  "comparison_group", "comparison_subgroup",
                  "comparability_confidence", "comparability_basis",
                  "comparable_for_matrix", "comparison_caveat",
                  "classification_source", "classification_source_table",
                  "classification_source_field",
                  "classification_source_value",
                  "classification_source_confidence")) {
      out[[col]] = .bundle_index(comparability, idx, col)
    }
  }
  for (col in c("metabolism_domain", "biosynthetic_family",
                "chemical_behavior", "comparison_scope",
                "comparison_group", "comparison_subgroup",
                "comparability_confidence", "comparability_basis",
                "comparable_for_matrix", "comparison_caveat")) {
    if (!col %in% names(out)) out[[col]] = NA_character_
  }

  first_cols = c(
    "species", "accepted_species_name", "genus", "family",
    "taxonomy_family_status", "matched_rank", "compound_id",
    "compound_name", "compound_name_clean", "CID", "pubchem_cid",
    "InChIKey", "SMILES", "MolecularFormula", "source_database",
    "source_record_id", "evidence_tier", "confidence", "occurrence_status",
    "occurrence_basis", "plant_part_group", "tissue_group", "method_group",
    "biological_context_status", "evidence_quality_score", "evidence_url",
    "doi", "pmid", "resolution_source", "identity_review_required",
    "context_known_record", "plant_part_known", "tissue_known",
    "method_known", "has_fingerprint", "metabolism_domain",
    "biosynthetic_family", "chemical_behavior", "comparison_scope",
    "comparison_group", "comparison_subgroup", "comparability_confidence",
    "comparability_basis", "comparable_for_matrix", "comparison_caveat",
    prop_cols, derived_cols
  )
  first_cols = unique(first_cols[first_cols %in% names(out)])
  out = out[, c(first_cols, setdiff(names(out), c(first_cols, "row_key"))),
            drop = FALSE]
  row.names(out) = NULL
  out
}

.bundle_membership_comparability = function(membership, derived) {
  derived_for_comp = derived
  if (is.data.frame(derived_for_comp) && nrow(derived_for_comp) > 0 &&
      "cid_key" %in% names(derived_for_comp)) {
    map = membership[!duplicated(.bundle_cid_key(membership)), ,
                     drop = FALSE]
    idx = match(derived_for_comp$cid_key, .bundle_cid_key(map))
    derived_for_comp$Query = .bundle_first_non_empty(
      .bundle_index(map, idx, "compound_name"),
      .bundle_index(map, idx, "compound_name_clean"),
      derived_for_comp$Query
    )
  }
  res = tryCatch(
    plantChemistryComparability(
      membership,
      categorate_result = list(DerivedGroups = derived_for_comp),
      min_confidence = "low"
    ),
    error = function(e) data.frame()
  )
  if (!is.data.frame(res)) return(data.frame())
  res
}

.bundle_enhance_pair_tanimoto = function(x) {
  x = as.data.frame(x, stringsAsFactors = FALSE)
  for (col in c("species_a", "species_b", "compound_pair_count")) {
    if (!col %in% names(x)) x[[col]] = NA
  }
  pair_a = pmin(x$species_a, x$species_b, na.rm = TRUE)
  pair_b = pmax(x$species_a, x$species_b, na.rm = TRUE)
  x$unordered_pair_key = paste(pair_a, pair_b, sep = " || ")
  count = suppressWarnings(as.numeric(x$compound_pair_count))
  x$support_tier = ifelse(is.na(count), "unknown",
                          ifelse(count >= 10000, "high",
                                 ifelse(count >= 1000, "moderate",
                                        ifelse(count >= 100, "low",
                                               "very_low"))))
  x$low_support_caution = .bundle_yes_no(!is.na(count) & count < 100)
  x$support_note = ifelse(
    is.na(count),
    "Compound-pair support was not reported.",
    ifelse(count < 100,
           "Interpret with caution: fewer than 100 compound-pair comparisons support this plant-pair summary.",
           "Pair summary has at least 100 compound-pair comparisons.")
  )
  dup = duplicated(x$unordered_pair_key)
  if (any(dup)) {
    x = x[!dup, , drop = FALSE]
  }
  row.names(x) = NULL
  x
}

.bundle_source_coverage_summary = function(coverage) {
  if (!is.data.frame(coverage) || nrow(coverage) < 1 ||
      !"Source" %in% names(coverage)) {
    return(data.frame(Source = character(), compounds_covered = integer(),
                      total_queries = integer(), total_records = numeric(),
                      coverage_fraction = numeric(),
                      interpretation_note = character(),
                      stringsAsFactors = FALSE))
  }
  if (!"Query" %in% names(coverage)) coverage$Query = NA_character_
  if (!"RecordCount" %in% names(coverage)) coverage$RecordCount = NA_integer_
  if (!"Present" %in% names(coverage)) {
    coverage$Present = suppressWarnings(as.numeric(coverage$RecordCount)) > 0
  }
  coverage$present_bool = .bundle_truthy(coverage$Present) |
    suppressWarnings(as.numeric(coverage$RecordCount)) > 0
  sources = sort(unique(.uaf_non_empty(coverage$Source)))
  rows = lapply(sources, function(src) {
    hit = coverage[coverage$Source == src, , drop = FALSE]
    total_queries = length(unique(.uaf_non_empty(hit$Query)))
    covered = length(unique(.uaf_non_empty(hit$Query[hit$present_bool])))
    records = sum(suppressWarnings(as.numeric(hit$RecordCount)), na.rm = TRUE)
    data.frame(
      Source = src,
      compounds_covered = covered,
      total_queries = total_queries,
      total_records = records,
      coverage_fraction = if (total_queries > 0) {
        round(covered / total_queries, 4)
      } else {
        NA_real_
      },
      interpretation_note = if (covered > 0) {
        "At least one record was returned for this source. Coverage does not imply complete annotation."
      } else {
        "No records were returned for this source in the exported enrichment batches."
      },
      stringsAsFactors = FALSE
    )
  })
  .bundle_bind(rows)
}

.bundle_species_summary = function(enriched) {
  cols = c("species", "genus", "family", "reported_compound_count",
           "fingerprinted_compound_count", "occurrence_record_count",
           "mean_evidence_quality_score", "plant_part_known_fraction",
           "context_known_record_fraction", "identity_review_compound_fraction",
           "dominant_comparison_scope", "dominant_comparison_group",
           "natural_product_fraction", "plant_occurring_fraction",
           "volatile_proxy_fraction", "lipophilic_fraction",
           "oxygenated_fraction", "nitrogenous_fraction",
           "sulfur_containing_fraction", "kegg_pathway_fraction",
           "safety_hazard_fraction", "literature_fraction",
           "source_database_count", "source_databases",
           "chemistry_data_quality_tier", "recommended_use_note")
  if (!is.data.frame(enriched) || nrow(enriched) < 1) {
    return(.bundle_empty(cols))
  }
  rows = lapply(split(enriched, enriched$species), function(x) {
    compound_ids = unique(.uaf_non_empty(x$compound_id))
    fingerprinted = unique(.uaf_non_empty(x$compound_id[
      .bundle_truthy(.bundle_col_or(x, "has_fingerprint", NA_character_))
    ]))
    occurrence_count = nrow(x)
    reported = length(compound_ids)
    fp_fraction = if (reported > 0) length(fingerprinted) / reported else 0
    data.frame(
      species = x$species[[1]],
      genus = .bundle_first_value(.bundle_col_or(x, "genus", NA_character_)),
      family = .bundle_first_value(.bundle_col_or(x, "family", NA_character_)),
      reported_compound_count = reported,
      fingerprinted_compound_count = length(fingerprinted),
      occurrence_record_count = occurrence_count,
      mean_evidence_quality_score = round(mean(suppressWarnings(as.numeric(
        .bundle_col_or(x, "evidence_quality_score", NA_real_)
      )), na.rm = TRUE), 4),
      plant_part_known_fraction = .bundle_fraction(
        .bundle_col_or(x, "plant_part_known", NA_character_)
      ),
      context_known_record_fraction = .bundle_fraction(
        .bundle_col_or(x, "context_known_record", NA_character_)
      ),
      identity_review_compound_fraction = .bundle_fraction(
        .bundle_col_or(x, "identity_review_required", NA_character_)
      ),
      dominant_comparison_scope = .bundle_mode(
        .bundle_col_or(x, "comparison_scope", NA_character_)
      ),
      dominant_comparison_group = .bundle_mode(
        .bundle_col_or(x, "comparison_group", NA_character_)
      ),
      natural_product_fraction = .bundle_fraction(
        .bundle_col_or(x, "is_natural_product", NA_character_)
      ),
      plant_occurring_fraction = .bundle_fraction(
        .bundle_col_or(x, "is_plant_occurring", NA_character_)
      ),
      volatile_proxy_fraction = .bundle_fraction(
        .bundle_col_or(x, "volatility_proxy", NA_character_),
        values = c("volatile_proxy", "volatile", "semivolatile")
      ),
      lipophilic_fraction = .bundle_fraction(
        .bundle_col_or(x, "lipophilicity_bin", NA_character_),
        values = c("lipophilic", "high_lipophilicity", "very_lipophilic")
      ),
      oxygenated_fraction = .bundle_fraction(
        .bundle_col_or(x, "oxygenated", NA_character_)
      ),
      nitrogenous_fraction = .bundle_fraction(
        .bundle_col_or(x, "nitrogenous", NA_character_)
      ),
      sulfur_containing_fraction = .bundle_fraction(
        .bundle_col_or(x, "sulfur_containing", NA_character_)
      ),
      kegg_pathway_fraction = .bundle_fraction(
        .bundle_col_or(x, "has_kegg_pathway", NA_character_)
      ),
      safety_hazard_fraction = .bundle_fraction(
        .bundle_col_or(x, "has_safety_hazard", NA_character_)
      ),
      literature_fraction = .bundle_fraction(
        .bundle_col_or(x, "has_literature", NA_character_)
      ),
      source_database_count = length(unique(.uaf_non_empty(
        .bundle_col_or(x, "source_database", NA_character_)
      ))),
      source_databases = .pubchem_collapse(unique(.uaf_non_empty(
        .bundle_col_or(x, "source_database", NA_character_)
      ))),
      chemistry_data_quality_tier = .bundle_species_quality_tier(
        reported, fp_fraction, occurrence_count
      ),
      recommended_use_note = .bundle_species_use_note(
        reported, fp_fraction
      ),
      stringsAsFactors = FALSE
    )
  })
  out = .bundle_bind(rows)
  out = out[, cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.bundle_missing_species = function(plant_list, enriched, metadata = NULL) {
  cols = c("species", "genus", "family", "matched_in_lotus",
           "reason_no_species_level_chemistry",
           "recommended_fallback_search", "recommended_manual_review",
           "notes")
  if (!is.data.frame(plant_list) || nrow(plant_list) < 1) {
    return(.bundle_empty(cols))
  }
  present = unique(.uaf_non_empty(enriched$species))
  missing = plant_list[!plant_list$species %in% present, , drop = FALSE]
  if (nrow(missing) < 1) return(.bundle_empty(cols))
  metadata = if (is.data.frame(metadata) && "species" %in% names(metadata)) {
    metadata[!duplicated(metadata$species), , drop = FALSE]
  } else {
    data.frame()
  }
  idx = if (nrow(metadata) > 0) match(missing$species, metadata$species) else
    rep(NA_integer_, nrow(missing))
  out = data.frame(
    species = missing$species,
    genus = .bundle_first_non_empty(.bundle_index(metadata, idx, "genus"),
                                    .bundle_genus(missing$species)),
    family = .bundle_index(metadata, idx, "family"),
    matched_in_lotus = "No",
    reason_no_species_level_chemistry =
      "No species-level public-source chemistry record was present in this bundle input.",
    recommended_fallback_search =
      "Review genus-level records, synonyms, KNApSAcK/NPASS/PubMed evidence, and project-specific literature before treating coverage as absent.",
    recommended_manual_review = "Yes",
    notes = ifelse(missing$species == "Elymus canadensis",
                   "Priority review: known-remediator benchmark missing species-level chemistry in this bundle.",
                   "Absence from this export is not biological absence of chemistry."),
    stringsAsFactors = FALSE
  )
  out[, cols, drop = FALSE]
}

.bundle_validation_overview = function(validation, batch_validation,
                                       validation_issues) {
  total_batches = if (is.data.frame(batch_validation)) nrow(batch_validation) else
    0L
  status = if (is.data.frame(batch_validation) &&
      "Status" %in% names(batch_validation)) {
    tolower(batch_validation$Status)
  } else {
    character()
  }
  issues = if (is.data.frame(validation_issues)) nrow(validation_issues) else 0L
  csv_status = if (is.list(validation) && is.data.frame(validation$Summary)) {
    validation$Summary$ExportReadyStatus[[1]]
  } else {
    "not_checked"
  }
  malformed = if (is.list(validation) &&
      is.data.frame(validation$CSVValidation) &&
      "BadRowCount" %in% names(validation$CSVValidation)) {
    sum(validation$CSVValidation$BadRowCount, na.rm = TRUE)
  } else {
    NA_integer_
  }
  blocking = identical(csv_status, "fail")
  data.frame(
    total_batches = total_batches,
    pass_batches = sum(status %in% c("pass", "ok"), na.rm = TRUE),
    warn_batches = sum(status %in% c("warn", "warning", "incomplete"),
                       na.rm = TRUE),
    error_batches = sum(status %in% c("error", "fail", "failed"),
                        na.rm = TRUE),
    total_warnings = issues,
    total_errors = sum(status %in% c("error", "fail", "failed"),
                       na.rm = TRUE),
    csv_parse_status = csv_status,
    malformed_csv_rows = malformed,
    export_ready_status = if (blocking) "fail" else "pass",
    validation_warning_classification =
      "Batch validation warnings should be reviewed. Duplicate analysis-key warnings are usually non-blocking when row-level evidence remains source-backed.",
    interpretation_note =
      "Export readiness only means the bundle is parseable and internally consistent; it does not mean public-source chemistry coverage is complete.",
    stringsAsFactors = FALSE
  )
}

.bundle_refresh_manifest = function(path, project_id = NULL) {
  existing = .bundle_read_manifest(path, required = FALSE)
  csv_files = sort(list.files(path, pattern = "[.]csv$", full.names = FALSE))
  csv_files = c("00_DataDictionary.csv",
                setdiff(csv_files, "00_DataDictionary.csv"))
  csv_files = csv_files[file.exists(file.path(path, csv_files))]
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  rows = lapply(csv_files, function(file_name) {
    table = .bundle_table_from_file(file_name, existing)
    dat = tryCatch(utils::read.csv(file.path(path, file_name),
                                   stringsAsFactors = FALSE,
                                   check.names = FALSE),
                   error = function(e) data.frame())
    data.frame(
      Table = table,
      SheetName = substr(gsub("[^A-Za-z0-9_]+", "_", table), 1, 31),
      FileName = file_name,
      RowCount = nrow(dat),
      ColumnCount = ncol(dat),
      Format = "csv",
      OutputPath = path,
      CreatedAt = created_at,
      ProjectID = .uaf_first_non_empty_text(project_id, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  manifest = .bundle_bind(rows)
  markdown = c("README.md", "METHODS_TEXT.md")
  md_rows = lapply(markdown[file.exists(file.path(path, markdown))],
                   function(file_name) {
    data.frame(
      Table = sub("[.].*$", "", file_name),
      SheetName = sub("[.].*$", "", file_name),
      FileName = file_name,
      RowCount = NA_integer_,
      ColumnCount = NA_integer_,
      Format = "markdown",
      OutputPath = path,
      CreatedAt = created_at,
      ProjectID = .uaf_first_non_empty_text(project_id, NA_character_),
      stringsAsFactors = FALSE
    )
  })
  manifest = .bundle_bind(c(list(manifest), md_rows))
  manifest$RowCount[manifest$Table == "ExportManifest"] = nrow(manifest)
  manifest$ColumnCount[manifest$Table == "ExportManifest"] = ncol(manifest)
  .categorate_write_csv_file(manifest,
                             file.path(path, "01_ExportManifest.csv"))
  row.names(manifest) = NULL
  manifest
}

.bundle_table_from_file = function(file_name, existing) {
  if (is.data.frame(existing) && nrow(existing) > 0 &&
      all(c("Table", "FileName") %in% names(existing))) {
    hit = match(file_name, existing$FileName)
    if (!is.na(hit)) return(existing$Table[[hit]])
  }
  known = c(
    "00_DataDictionary.csv" = "DataDictionary",
    "01_ExportManifest.csv" = "ExportManifest",
    "02_BatchSummary.csv" = "BatchSummary",
    "03_PlantCompoundMembership.csv" = "PlantCompoundMembership",
    "03b_PlantCompoundMembershipEnriched.csv" =
      "PlantCompoundMembershipEnriched",
    "04_PlantPairTanimotoSummary.csv" = "PlantPairTanimotoSummary",
    "05_ResolvedCompounds.csv" = "ResolvedCompounds",
    "06_PubChemFingerprints.csv" = "PubChemFingerprints",
    "07_FileReferences.csv" = "FileReferences",
    "08_ChemicalTraitSummary.csv" = "ChemicalTraitSummary",
    "09_DerivedGroups.csv" = "DerivedGroups",
    "10_PubChemProperties.csv" = "PubChemProperties",
    "11_SourceCoverage.csv" = "SourceCoverage",
    "11b_SourceCoverageSummary.csv" = "SourceCoverageSummary",
    "12_ValidationSummary.csv" = "ValidationSummary",
    "12b_ValidationOverview.csv" = "ValidationOverview",
    "13_ValidationIssues.csv" = "ValidationIssues",
    "14_PlantChemistrySummary.csv" = "PlantChemistrySummary",
    "15_PlantChemistryMissingSpecies.csv" = "PlantChemistryMissingSpecies"
  )
  if (file_name %in% names(known)) return(unname(known[[file_name]]))
  sub("[.]csv$", "", sub("^[0-9]+[A-Za-z]?_", "", file_name))
}

.bundle_data_dictionary = function(path, manifest) {
  csv_manifest = manifest[manifest$Format == "csv", , drop = FALSE]
  rows = list()
  for (i in seq_len(nrow(csv_manifest))) {
    file = file.path(path, csv_manifest$FileName[[i]])
    if (!file.exists(file)) next
    dat = tryCatch(utils::read.csv(file, nrows = 50,
                                   stringsAsFactors = FALSE,
                                   check.names = FALSE),
                   error = function(e) data.frame())
    if (ncol(dat) < 1) next
    for (col in names(dat)) {
      rows[[length(rows) + 1L]] = data.frame(
        table_name = csv_manifest$Table[[i]],
        file_name = csv_manifest$FileName[[i]],
        table_purpose = .bundle_table_purpose(csv_manifest$Table[[i]]),
        column_name = col,
        data_type = class(dat[[col]])[[1]],
        required = .bundle_yes_no(.bundle_column_required(
          csv_manifest$Table[[i]], col
        )),
        allowed_values = .bundle_allowed_values(col),
        interpretation_notes = .bundle_column_note(col),
        known_limitations = .bundle_column_limitation(col),
        stringsAsFactors = FALSE
      )
    }
  }
  .bundle_bind(rows)
}

.bundle_write_readme = function(path, project_id = NULL) {
  lines = c(
    "# Plant Chemistry Analysis Bundle",
    "",
    paste0("Project: ", .uaf_first_non_empty_text(project_id, "unspecified")),
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
    "",
    "This bundle is a canonical uafR handoff for plant chemistry analysis. It combines species-level public-source plant-compound memberships, PubChem fingerprint resolution, plant-pair Tanimoto summaries, categorate annotation tables, source coverage, validation outputs, and analysis-ready summaries.",
    "",
    "Primary analysis files:",
    "- `03_PlantCompoundMembership.csv`: source-backed plant-compound occurrence records.",
    "- `03b_PlantCompoundMembershipEnriched.csv`: membership joined to identity, fingerprint, PubChem property, derived group, context, and comparability fields.",
    "- `04_PlantPairTanimotoSummary.csv`: plant-pair structural similarity summaries with unordered pair keys and support tiers.",
    "- `14_PlantChemistrySummary.csv`: species-level chemistry coverage and analysis-readiness summary.",
    "- `15_PlantChemistryMissingSpecies.csv`: supplied species with no species-level chemistry in this bundle input.",
    "",
    "QA and provenance files:",
    "- `00_DataDictionary.csv`: table and column definitions.",
    "- `01_ExportManifest.csv`: file names, row counts, column counts, and creation metadata.",
    "- `11b_SourceCoverageSummary.csv`: source-level coverage summary.",
    "- `12b_ValidationOverview.csv`: export-level validation status.",
    "- `13_ValidationIssues.csv`: row-level categorate validation issues.",
    "",
    "Interpretation limits:",
    "- These records are public-source, species-level reported chemistry records, not direct measurements from project samples.",
    "- Absence of a compound record is not evidence that a plant lacks the compound.",
    "- Most public records lack plant-part context; missing context should not be interpreted as whole-plant evidence.",
    "- LOTUS or other database occurrence is not proof of pathway activity, ecological function, or phytoremediation performance.",
    "- PubChem structural similarity is not functional equivalence.",
    "- Categorate and PubChem annotations are source evidence for interpretation, not experimental confirmation.",
    "- Phytoremediation claims require contaminant-specific literature or experimental evidence."
  )
  writeLines(lines, file.path(path, "README.md"), useBytes = TRUE)
}

.bundle_write_methods = function(path, project_id = NULL) {
  lines = c(
    "# Methods And Provenance Text",
    "",
    paste0("Project: ", .uaf_first_non_empty_text(project_id, "unspecified")),
    "",
    "Plant-compound memberships were exported from cached uafR plant chemistry outputs. Occurrence rows should be interpreted as public-source species-level chemistry evidence. Evidence tiers separate direct species database records from weaker candidate or fallback evidence when those fields are present.",
    "",
    "Compound identities were carried forward from source-backed uafR resolution fields. PubChem CIDs, fingerprints, structural strings, molecular formulae, and computed/annotated properties were joined only when present in cached PubChem or categorate outputs. No SMILES, InChIKeys, formulas, CIDs, or abundance-like values were fabricated.",
    "",
    "Plant-pair Tanimoto summaries were computed from PubChem Fingerprint2D compound fingerprints through the uafR plant chemistry Tanimoto workflow. The summary reports pair-level aggregate structural similarity and support counts. Structural similarity does not imply equivalent biological role, contaminant response, or pathway activity.",
    "",
    "Comparability fields classify compounds into conservative chemistry scopes using source-backed uafR derived groups, chemical class fields, and transparent fallback rules. These fields are intended to prevent analyses from mixing unlike chemistry, such as primary metabolites, specialized metabolites, volatile-specialized chemistry, lipids/fatty acids, hormone signals, and unknown chemistry.",
    "",
    "Validation outputs report parseability, manifest consistency, required-column completeness, source coverage, and categorate warning counts. A parseable, internally consistent bundle is suitable for downstream analysis, but manuscript claims must still be limited to the strength of the public-source evidence."
  )
  writeLines(lines, file.path(path, "METHODS_TEXT.md"), useBytes = TRUE)
}

.bundle_validate_csv_file = function(file, manifest, use_python = FALSE,
                                     use_pandas = FALSE) {
  file_name = basename(file)
  expected_rows = expected_cols = NA_integer_
  if (is.data.frame(manifest) && nrow(manifest) > 0 &&
      "FileName" %in% names(manifest)) {
    hit = match(file_name, manifest$FileName)
    if (!is.na(hit)) {
      expected_rows = suppressWarnings(as.integer(manifest$RowCount[[hit]]))
      expected_cols = suppressWarnings(as.integer(manifest$ColumnCount[[hit]]))
    }
  }
  warnings = character()
  dat = withCallingHandlers(
    tryCatch(utils::read.csv(file, stringsAsFactors = FALSE,
                             check.names = FALSE),
             error = function(e) e),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  r_ok = is.data.frame(dat)
  actual_rows = if (r_ok) nrow(dat) else NA_integer_
  actual_cols = if (r_ok) ncol(dat) else NA_integer_
  fields = tryCatch(utils::count.fields(file, sep = ",", quote = "\"",
                                        blank.lines.skip = FALSE),
                    error = function(e) NA_integer_)
  header_cols = if (length(fields) > 0) fields[[1]] else NA_integer_
  bad = if (length(fields) > 0 && !all(is.na(fields))) {
    which(fields != header_cols)
  } else {
    integer()
  }
  py = if (isTRUE(use_python) || isTRUE(use_pandas)) {
    .bundle_python_csv_check(file, expected_cols, use_pandas)
  } else {
    list(csv_status = "not_checked", pandas_status = "not_checked",
         message = NA_character_)
  }
  extra_index = r_ok && ncol(dat) > 0 &&
    names(dat)[[1]] %in% c("X", "...1", "row.names")
  status = if (!r_ok || length(bad) > 0 ||
      (!is.na(expected_rows) && !identical(actual_rows, expected_rows)) ||
      (!is.na(expected_cols) && !identical(actual_cols, expected_cols)) ||
      identical(py$csv_status, "fail") || identical(py$pandas_status, "fail")) {
    "fail"
  } else if (length(warnings) > 0 || extra_index) {
    "warn"
  } else {
    "pass"
  }
  data.frame(
    FileName = file_name,
    ExpectedRows = expected_rows,
    ActualRows = actual_rows,
    ExpectedColumns = expected_cols,
    ActualColumns = actual_cols,
    HeaderColumns = header_cols,
    BadRowCount = length(bad),
    BadRows = .pubchem_collapse(utils::head(bad, 20)),
    RReadStatus = if (r_ok) "pass" else "fail",
    PythonCsvStatus = py$csv_status,
    PandasStatus = py$pandas_status,
    ExtraIndexColumn = .bundle_yes_no(extra_index),
    WarningMessages = .pubchem_collapse(unique(warnings)),
    ParserMessage = py$message,
    Status = status,
    stringsAsFactors = FALSE
  )
}

.bundle_python_csv_check = function(file, expected_cols, use_pandas) {
  python = Sys.which("python3")
  if (!nzchar(python)) {
    return(list(csv_status = "not_available",
                pandas_status = if (isTRUE(use_pandas)) "not_available" else
                  "not_checked",
                message = "python3 was not available."))
  }
  script = tempfile(fileext = ".py")
  writeLines(c(
    "import csv, json, sys",
    "path = sys.argv[1]",
    "expected = None if sys.argv[2] == 'NA' else int(sys.argv[2])",
    "want_pandas = sys.argv[3] == 'TRUE'",
    "out = {'csv_status':'pass','pandas_status':'not_checked','message':''}",
    "try:",
    "    with open(path, newline='', encoding='utf-8') as f:",
    "        reader = csv.reader(f)",
    "        header = next(reader)",
    "        n = len(header)",
    "        bad = []",
    "        for i, row in enumerate(reader, start=2):",
    "            if len(row) != n:",
    "                bad.append(i)",
    "        if expected is not None and n != expected:",
    "            bad.append(1)",
    "        if bad:",
    "            out['csv_status'] = 'fail'",
    "            out['message'] = 'bad csv rows: ' + ','.join(map(str, bad[:20]))",
    "except Exception as e:",
    "    out['csv_status'] = 'fail'",
    "    out['message'] = str(e)",
    "if want_pandas:",
    "    try:",
    "        import pandas as pd",
    "        pd.read_csv(path)",
    "        out['pandas_status'] = 'pass'",
    "    except Exception as e:",
    "        out['pandas_status'] = 'fail'",
    "        out['message'] = (out.get('message','') + ' pandas: ' + str(e)).strip()",
    "print(json.dumps(out))"
  ), script, useBytes = TRUE)
  args = c(script, file, ifelse(is.na(expected_cols), "NA",
                                as.character(expected_cols)),
           ifelse(isTRUE(use_pandas), "TRUE", "FALSE"))
  res = tryCatch(system2(python, args = args, stdout = TRUE, stderr = TRUE),
                 error = function(e) character())
  out = tryCatch(jsonlite::fromJSON(paste(res, collapse = "\n")),
                 error = function(e) list(csv_status = "fail",
                                          pandas_status = "fail",
                                          message = paste(res, collapse = " ")))
  unlink(script)
  out
}

.bundle_validate_required_columns = function(path, manifest) {
  required = list(
    PlantCompoundMembership = c("species", "compound_id", "compound_name",
                                "source_database", "evidence_tier"),
    PlantCompoundMembershipEnriched = c(
      "species", "compound_id", "compound_name", "comparison_scope",
      "comparison_group", "has_fingerprint"
    ),
    PlantPairTanimotoSummary = c("species_a", "species_b",
                                 "compound_pair_count", "mean_tanimoto",
                                 "unordered_pair_key", "support_tier"),
    PlantChemistrySummary = c("species", "reported_compound_count",
                              "fingerprinted_compound_count",
                              "chemistry_data_quality_tier")
  )
  rows = lapply(names(required), function(table) {
    if (!is.data.frame(manifest) || !"Table" %in% names(manifest) ||
        !table %in% manifest$Table) {
      return(data.frame(
        Table = table,
        RequiredColumns = .pubchem_collapse(required[[table]]),
        MissingColumns = NA_character_,
        Status = "not_present",
        stringsAsFactors = FALSE
      ))
    }
    dat = .bundle_read_table(path, manifest, table)
    missing = setdiff(required[[table]], names(dat))
    data.frame(
      Table = table,
      RequiredColumns = .pubchem_collapse(required[[table]]),
      MissingColumns = .pubchem_collapse(missing),
      Status = if (length(missing) > 0) "fail" else "pass",
      stringsAsFactors = FALSE
    )
  })
  .bundle_bind(rows)
}

.bundle_prepare_cid_table = function(x) {
  if (!is.data.frame(x) || nrow(x) < 1) {
    return(data.frame(cid_key = character(), stringsAsFactors = FALSE))
  }
  x = as.data.frame(x, stringsAsFactors = FALSE)
  x$cid_key = .bundle_cid_key(x)
  x = x[!is.na(x$cid_key) & x$cid_key != "", , drop = FALSE]
  x = x[!duplicated(x$cid_key), , drop = FALSE]
  row.names(x) = NULL
  x
}

.bundle_cid_key = function(x) {
  if (!is.data.frame(x) || nrow(x) < 1) return(character())
  value = rep(NA_character_, nrow(x))
  for (col in c("pubchem_cid", "CID", "source_CID")) {
    if (col %in% names(x)) {
      value = .bundle_first_non_empty(value, as.character(x[[col]]))
    }
  }
  if ("Query" %in% names(x)) {
    query_cid = sub("^cid:", "", as.character(x$Query), ignore.case = TRUE)
    query_cid[!grepl("^cid:", as.character(x$Query), ignore.case = TRUE)] =
      NA_character_
    value = .bundle_first_non_empty(value, query_cid)
  }
  value = gsub("[^0-9]", "", value)
  value[value == ""] = NA_character_
  value
}

.bundle_unique_by = function(x, by) {
  if (!is.data.frame(x) || nrow(x) < 1 || !by %in% names(x)) {
    return(data.frame())
  }
  x = as.data.frame(x, stringsAsFactors = FALSE)
  x = x[!is.na(x[[by]]) & x[[by]] != "", , drop = FALSE]
  x = x[!duplicated(x[[by]]), , drop = FALSE]
  row.names(x) = NULL
  x
}

.bundle_index = function(x, idx, col) {
  if (!is.data.frame(x) || nrow(x) < 1 || !col %in% names(x)) {
    return(rep(NA_character_, length(idx)))
  }
  out = rep(NA_character_, length(idx))
  ok = !is.na(idx) & idx >= 1 & idx <= nrow(x)
  out[ok] = as.character(x[[col]][idx[ok]])
  out
}

.bundle_col_or = function(x, col, default) {
  if (is.data.frame(x) && col %in% names(x)) return(x[[col]])
  rep(default, if (is.data.frame(x)) nrow(x) else 0L)
}

.bundle_first_non_empty = function(...) {
  values = list(...)
  if (length(values) < 1) return(character())
  n = max(vapply(values, length, integer(1)))
  values = lapply(values, function(x) {
    x = as.character(x)
    if (length(x) == 1 && n > 1) x = rep(x, n)
    if (length(x) < n) x = c(x, rep(NA_character_, n - length(x)))
    x
  })
  out = rep(NA_character_, n)
  for (value in values) {
    hit = (is.na(out) | out == "") & !is.na(value) & value != ""
    out[hit] = value[hit]
  }
  out
}

.bundle_first_value = function(x) {
  hit = .uaf_non_empty(x)
  if (length(hit) < 1) return(NA_character_)
  hit[[1]]
}

.bundle_occurrence_key = function(x) {
  paste(.bundle_col_or(x, "species", NA_character_),
        .bundle_col_or(x, "compound_name_clean", NA_character_),
        .bundle_col_or(x, "source_database", NA_character_),
        .bundle_col_or(x, "source_record_id", NA_character_),
        sep = "\r")
}

.bundle_known = function(x) {
  x = tolower(.bundle_squish(as.character(x)))
  !is.na(x) & x != "" & !x %in% c("unknown", "unspecified", "not_reported",
                                  "not reported", "na", "none")
}

.bundle_truthy = function(x) {
  x = tolower(.bundle_squish(as.character(x)))
  x %in% c("true", "t", "yes", "y", "1", "present", "high", "medium")
}

.bundle_yes_no = function(x) ifelse(isTRUE(x) | (!is.na(x) & x), "Yes", "No")

.bundle_fraction = function(x, values = NULL) {
  if (length(x) < 1) return(NA_real_)
  if (is.null(values)) {
    val = .bundle_truthy(x)
  } else {
    val = tolower(.bundle_squish(as.character(x))) %in% tolower(values)
  }
  round(mean(val, na.rm = TRUE), 4)
}

.bundle_mode = function(x) {
  x = .uaf_non_empty(x)
  x = x[!tolower(x) %in% c("unknown", "none")]
  if (length(x) < 1) return(NA_character_)
  tab = sort(table(x), decreasing = TRUE)
  names(tab)[[1]]
}

.bundle_species_quality_tier = function(reported, fp_fraction,
                                        occurrence_count) {
  if (reported < 1) return("no_species_level_public_profile")
  if (reported >= 50 && fp_fraction >= 0.8 && occurrence_count >= 50) {
    return("robust_public_species_profile")
  }
  if (reported >= 10 && fp_fraction >= 0.5) {
    return("moderate_public_species_profile")
  }
  "limited_public_species_profile"
}

.bundle_species_use_note = function(reported, fp_fraction) {
  if (reported < 1) {
    return("No species-level public chemistry records were included; use only for missing-coverage reporting.")
  }
  if (fp_fraction < 0.5) {
    return("Use cautiously in structure-based analyses because fewer than half of reported compounds have fingerprints.")
  }
  "Suitable for source-backed exploratory chemistry comparison with stated public-data limitations."
}

.bundle_table_purpose = function(table) {
  switch(table,
         DataDictionary = "Defines exported tables and columns.",
         ExportManifest = "Lists bundle files, row counts, and provenance.",
         PlantCompoundMembership = "Raw source-backed plant-compound occurrence records.",
         PlantCompoundMembershipEnriched = "Analysis-ready occurrence records joined to identity, property, context, and comparability fields.",
         PlantPairTanimotoSummary = "Plant-pair structural similarity summaries from compound fingerprints.",
         PlantChemistrySummary = "Species-level chemistry coverage and analysis-readiness summary.",
         PlantChemistryMissingSpecies = "Supplied species with no species-level chemistry record in this bundle.",
         SourceCoverageSummary = "Source-level coverage summary from categorate batches.",
         ValidationOverview = "Bundle-level validation and export readiness summary.",
         "uafR plant chemistry bundle table.")
}

.bundle_column_required = function(table, col) {
  req = list(
    PlantCompoundMembership = c("species", "compound_id", "compound_name"),
    PlantCompoundMembershipEnriched = c("species", "compound_id",
                                        "compound_name", "comparison_scope"),
    PlantPairTanimotoSummary = c("species_a", "species_b",
                                 "unordered_pair_key"),
    PlantChemistrySummary = c("species", "reported_compound_count")
  )
  table %in% names(req) && col %in% req[[table]]
}

.bundle_allowed_values = function(col) {
  switch(col,
         comparable_for_matrix = "Yes; No",
         has_fingerprint = "Yes; No",
         context_known_record = "Yes; No",
         plant_part_known = "Yes; No",
         tissue_known = "Yes; No",
         method_known = "Yes; No",
         matched_in_lotus = "Yes; No",
         support_tier = "very_low; low; moderate; high; unknown",
         ExportReadyStatus = "pass; warn; fail",
         "")
}

.bundle_column_note = function(col) {
  switch(col,
         comparison_scope = "Use to keep downstream comparisons biologically comparable.",
         comparison_group = "Higher-level chemistry grouping within a comparison scope.",
         evidence_tier = "Source evidence level; literature co-mentions should not be treated as confirmed occurrence unless curated.",
         plant_part_group = "Often missing in public sources; absence means not reported in this export.",
         mean_tanimoto = "Mean PubChem Fingerprint2D Tanimoto similarity across supported compound pairs.",
         support_tier = "Support tier is based on compound_pair_count and guides caution.",
         chemistry_data_quality_tier = "Coverage tier for exploratory analysis, not a claim of biological completeness.",
         "")
}

.bundle_column_limitation = function(col) {
  switch(col,
         is_plant_occurring = "Source annotation presence does not prove occurrence in the project sample.",
         has_kegg_pathway = "Pathway context is not evidence of pathway activity in the plant.",
         has_safety_hazard = "Hazard annotation is screening evidence, not risk assessment.",
         mean_tanimoto = "Structural similarity is not functional equivalence.",
         family = "Only populated when supplied by source data or metadata; not inferred by uafR during bundle finalization.",
         "")
}

.bundle_empty = function(cols) {
  out = as.data.frame(stats::setNames(rep(list(character()), length(cols)),
                                      cols),
                      stringsAsFactors = FALSE)
  out
}

.bundle_bind = function(rows) {
  rows = rows[vapply(rows, is.data.frame, logical(1))]
  if (length(rows) < 1) return(data.frame())
  cols = unique(unlist(lapply(rows, names), use.names = FALSE))
  aligned = lapply(rows, function(x) {
    for (col in setdiff(cols, names(x))) x[[col]] = NA
    x[, cols, drop = FALSE]
  })
  out = do.call(rbind, aligned)
  row.names(out) = NULL
  out
}

.bundle_genus = function(species) {
  species = .bundle_squish(species)
  vapply(strsplit(species, "[[:space:]]+"), function(x) {
    if (length(x) < 1) return(NA_character_)
    x[[1]]
  }, character(1))
}

.bundle_squish = function(x) {
  x = as.character(x)
  x = gsub("[[:space:]]+", " ", x)
  trimws(x)
}
