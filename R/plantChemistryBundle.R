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
#' @param plant_compound_pair_tanimoto Optional plant-compound pair Tanimoto
#' data frame or CSV/CSV.GZ path used to build comparable scope/group summaries
#' when `include_comparable_tanimoto = TRUE`.
#' @param project_id Optional project label written to documentation.
#' @param overwrite Logical. If `FALSE`, existing finalization outputs are not
#' overwritten.
#' @param include_feature_exports Logical. If `TRUE`, write model-ready species
#' feature matrices and species quality metadata.
#' @param include_comparable_tanimoto Logical. If `TRUE`, write Tanimoto
#' summaries filtered to matched comparable chemistry scopes and groups.
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
                                                plant_compound_pair_tanimoto = NULL,
                                                project_id = NULL,
                                                overwrite = TRUE,
                                                include_feature_exports = TRUE,
                                                include_comparable_tanimoto = FALSE,
                                                validate_export = TRUE,
                                                max_cell_chars = 30000) {
  path = .uaf_non_empty(path)
  if (length(path) != 1 || !dir.exists(path)) {
    stop("`path` must be one existing CSV bundle directory.", call. = FALSE)
  }
  path = normalizePath(path, winslash = "/", mustWork = FALSE)
  manifest = .bundle_read_manifest(path)
  metadata = .bundle_optional_frame(metadata, "metadata")
  if (is.data.frame(metadata)) metadata = standardizePlantMetadata(metadata)
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

  evidence_grades = plantOccurrenceEvidenceGrade(enriched)
  .bundle_write_table(path, "16_EvidenceGradeSummary.csv",
                      .bundle_evidence_grade_summary(evidence_grades),
                      overwrite = overwrite, max_cell_chars = max_cell_chars)
  .bundle_write_table(path, "17_ReviewRequiredOccurrences.csv",
                      .bundle_review_required_occurrences(enriched,
                                                          evidence_grades),
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  if (is.data.frame(pair_tanimoto) && nrow(pair_tanimoto) > 0) {
    pair_tanimoto = .bundle_enhance_pair_tanimoto(pair_tanimoto)
    .bundle_write_table(path, "04_PlantPairTanimotoSummary.csv",
                        pair_tanimoto, overwrite = TRUE,
                        max_cell_chars = max_cell_chars)
  }
  if (isTRUE(include_comparable_tanimoto)) {
    comparable = plantComparableTanimotoSummary(
      plant_pair_tanimoto = pair_tanimoto,
      plant_compound_pair_tanimoto = plant_compound_pair_tanimoto,
      membership = enriched
    )
    .bundle_write_table(path, "22_ComparableScopeTanimotoSummary.csv",
                        comparable$ScopeFiltered, overwrite = overwrite,
                        max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "23_ComparableGroupTanimotoSummary.csv",
                        comparable$GroupFiltered, overwrite = overwrite,
                        max_cell_chars = max_cell_chars)
  }

  source_summary = .bundle_source_coverage_summary(coverage)
  .bundle_write_table(path, "11b_SourceCoverageSummary.csv", source_summary,
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  species_summary = .bundle_species_summary(enriched)
  .bundle_write_table(path, "14_PlantChemistrySummary.csv", species_summary,
                      overwrite = overwrite, max_cell_chars = max_cell_chars)

  if (isTRUE(include_feature_exports)) {
    feature_universe = if (is.data.frame(plant_list) &&
                           nrow(plant_list) > 0) {
      plant_list$species
    } else {
      enriched$species
    }
    features = exportPlantChemistryFeatureSet(
      enriched,
      species_universe = feature_universe,
      plant_metadata = metadata
    )
    .bundle_write_table(path, "18_FeatureComparisonGroupCountMatrix.csv",
                        features$ComparisonGroupCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "19_FeatureComparisonScopeCountMatrix.csv",
                        features$ComparisonScopeCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "20_FeatureSourceCoverageMatrix.csv",
                        features$SourceCoverageMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "21_FeatureSpeciesMetadata.csv",
                        features$SpeciesMetadata,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "24_FeatureEvidenceGradeCountMatrix.csv",
                        features$EvidenceGradeCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "25_FeaturePlantPartCountMatrix.csv",
                        features$PlantPartCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "26_FeatureTissueCountMatrix.csv",
                        features$TissueCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    .bundle_write_table(path, "27_FeatureMethodCountMatrix.csv",
                        features$MethodCountMatrix,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
    feature_manifest = .bundle_finalize_feature_manifest(features$Manifest)
    .bundle_write_table(path, "28_FeatureMatrixManifest.csv",
                        feature_manifest,
                        overwrite = overwrite, max_cell_chars = max_cell_chars)
  }

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
#' consistency, manifest row/column agreement, portable artifact paths,
#' file sizes and checksums, accidental row-index columns, and required columns
#' in the main analysis tables. Python `csv.reader` and pandas checks are
#' optional so package tests do not depend on a Python installation, but the
#' function records whether those checks were run.
#'
#' @param path CSV bundle directory.
#' @param use_python Logical. If `TRUE`, also validate every CSV with Python's
#' standard `csv.reader` when `python3` is available.
#' @param use_pandas Logical. If `TRUE`, also validate every CSV with
#' `pandas.read_csv()` when Python and pandas are available.
#'
#' @return A list with `Summary`, `CSVValidation`, `RequiredColumns`,
#' `ArtifactValidation`, `ManifestReferences`, and `Manifest` tables.
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
  artifacts = .bundle_validate_export_artifacts(path, manifest)
  manifest_references = .bundle_validate_manifest_references(path, manifest)
  blocking = any(csv_validation$Status == "fail") ||
    any(required$Status == "fail") ||
    any(artifacts$Status == "fail") ||
    any(manifest_references$Status == "fail")
  warnings = any(csv_validation$Status == "warn") ||
    any(required$Status == "warn") ||
    any(artifacts$Status == "warn") ||
    any(manifest_references$Status == "warn")
  summary = data.frame(
    CSVFileCount = length(csv_files),
    CSVPassCount = sum(csv_validation$Status == "pass", na.rm = TRUE),
    CSVWarnCount = sum(csv_validation$Status == "warn", na.rm = TRUE),
    CSVFailCount = sum(csv_validation$Status == "fail", na.rm = TRUE),
    RequiredColumnFailCount = sum(required$Status == "fail", na.rm = TRUE),
    ArtifactWarnCount = sum(artifacts$Status == "warn", na.rm = TRUE),
    ArtifactFailCount = sum(artifacts$Status == "fail", na.rm = TRUE),
    ManifestReferenceFailCount = sum(manifest_references$Status == "fail",
                                     na.rm = TRUE),
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
    ArtifactValidation = artifacts,
    ManifestReferences = manifest_references,
    Manifest = manifest
  )
}

#' Standardize plant metadata for chemistry projects
#'
#' @description
#' `standardizePlantMetadata()` normalizes user-supplied plant metadata for
#' joins into plant chemistry bundles. It derives genus from the supplied
#' species name when genus is absent, but it does not infer family or accepted
#' names. Missing family and accepted-name coverage are reported through status
#' columns so downstream bundles stay transparent.
#'
#' @param x Data frame, CSV path, or character vector of species names.
#' @param species_col Optional species column name.
#' @param accepted_species_col Optional accepted species column name.
#' @param genus_col Optional genus column name.
#' @param family_col Optional family column name.
#' @param group_col Optional project group column name.
#' @param role_col Optional project role/label column name.
#' @param metadata_source Source label written to `metadata_source`.
#'
#' @return Standardized plant metadata data frame.
#'
#' @export
standardizePlantMetadata = function(x,
                                    species_col = NULL,
                                    accepted_species_col = NULL,
                                    genus_col = NULL,
                                    family_col = NULL,
                                    group_col = NULL,
                                    role_col = NULL,
                                    metadata_source = "user_supplied") {
  if (is.character(x) && length(x) == 1 && file.exists(x)) {
    x = utils::read.csv(x, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (is.character(x)) {
    x = data.frame(species = x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) {
    stop("`x` must be a metadata data frame, CSV path, or species vector.",
         call. = FALSE)
  }
  if (nrow(x) < 1) {
    return(data.frame(
      species = character(), accepted_species_name = character(),
      genus = character(), family = character(), group = character(),
      role = character(), metadata_source = character(),
      accepted_name_status = character(), family_status = character(),
      stringsAsFactors = FALSE
    ))
  }
  raw_names = names(x)
  names(x) = .plant_normalize_column_names(names(x))
  species_col = .metadata_detect_col(x, species_col,
                                     c("species", "scientific_name",
                                       "taxon", "plant", "plant_species"))
  if (is.null(species_col)) {
    stop("Plant metadata requires a species/scientific-name column.",
         call. = FALSE)
  }
  accepted_species_col = .metadata_detect_col(
    x, accepted_species_col,
    c("accepted_species_name", "accepted_name", "accepted_taxon",
      "accepted_scientific_name")
  )
  genus_col = .metadata_detect_col(x, genus_col, c("genus"))
  family_col = .metadata_detect_col(x, family_col, c("family", "plant_family"))
  group_col = .metadata_detect_col(x, group_col, c("group", "class", "label"))
  role_col = .metadata_detect_col(x, role_col, c("role", "function",
                                                "project_role"))

  species = .bundle_squish(x[[species_col]])
  accepted = if (!is.null(accepted_species_col)) {
    .bundle_squish(x[[accepted_species_col]])
  } else {
    rep(NA_character_, nrow(x))
  }
  genus = if (!is.null(genus_col)) .bundle_squish(x[[genus_col]]) else
    .bundle_genus(species)
  genus = .bundle_first_non_empty(genus, .bundle_genus(species))
  family = if (!is.null(family_col)) .bundle_squish(x[[family_col]]) else
    rep(NA_character_, nrow(x))
  group = if (!is.null(group_col)) .bundle_squish(x[[group_col]]) else
    rep(NA_character_, nrow(x))
  role = if (!is.null(role_col)) .bundle_squish(x[[role_col]]) else
    rep(NA_character_, nrow(x))

  out = data.frame(
    species = species,
    accepted_species_name = accepted,
    genus = genus,
    family = family,
    group = group,
    role = role,
    metadata_source = metadata_source,
    accepted_name_status = ifelse(.bundle_known(accepted),
                                  "accepted_name_supplied",
                                  "accepted_name_not_supplied"),
    family_status = ifelse(.bundle_known(family), "family_supplied",
                           "family_not_supplied"),
    stringsAsFactors = FALSE
  )
  keep = !is.na(out$species) & out$species != ""
  out = out[keep, , drop = FALSE]
  out = out[!duplicated(out$species), , drop = FALSE]
  row.names(out) = NULL
  source_cols = c(species = species_col,
                  accepted_species_name = accepted_species_col,
                  genus = genus_col,
                  family = family_col,
                  group = group_col,
                  role = role_col)
  source_cols[vapply(source_cols, is.null, logical(1))] = NA_character_
  attr(out, "source_columns") = unlist(source_cols, use.names = TRUE)
  attr(out, "raw_column_names") = raw_names
  out
}

#' Grade plant-compound occurrence evidence
#'
#' @description
#' `plantOccurrenceEvidenceGrade()` turns source evidence fields into a
#' transparent analysis tier. Grades are conservative and do not claim that a
#' compound was measured in project samples.
#'
#' @param x Plant phytochemistry result, occurrence table, or enriched
#' plant-compound membership table.
#' @param review_table Optional review table created by uafR review workflows.
#'
#' @return Row-level evidence grade table.
#'
#' @export
plantOccurrenceEvidenceGrade = function(x, review_table = NULL) {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  if (!is.data.frame(occurrences) || nrow(occurrences) < 1) {
    return(.bundle_empty(.evidence_grade_cols()))
  }
  occurrences = as.data.frame(occurrences, stringsAsFactors = FALSE)
  for (col in c("species", "compound_id", "compound_name",
                "compound_name_clean", "source_database", "source_record_id",
                "evidence_tier", "matched_rank", "occurrence_status",
                "confidence", "analysis_ready", "identity_review_required",
                "SMILES", "InChIKey", "CID", "pubchem_cid",
                "comparable_for_matrix", "comparison_scope",
                "comparison_group", "curation_flag", "evidence_url",
                "pmid", "doi")) {
    if (!col %in% names(occurrences)) occurrences[[col]] = NA_character_
  }
  review = .evidence_review_decisions(occurrences, review_table)
  rows = lapply(seq_len(nrow(occurrences)), function(i) {
    row = occurrences[i, , drop = FALSE]
    decision = review$review_decision[[i]]
    grade = .evidence_grade_for_row(row, decision)
    structure_resolved = .bundle_known(row$SMILES) |
      .bundle_known(row$InChIKey) | .bundle_known(row$CID) |
      .bundle_known(row$pubchem_cid)
    source_backed = !grade$grade %in%
      c("pubtator_pubmed_candidate_only", "unresolved_or_review_required",
        "excluded_by_review")
    review_required = grade$grade %in%
      c("pubtator_pubmed_candidate_only", "unresolved_or_review_required",
        "source_backed_genus_family_fallback") |
      .bundle_truthy(row$identity_review_required) |
      grepl("review", .bundle_squish(row$curation_flag), ignore.case = TRUE)
    comparable = .bundle_truthy(row$comparable_for_matrix) &&
      !tolower(.bundle_squish(row$comparison_scope)) %in% c("", "unknown")
    data.frame(
      row_id = i,
      species = row$species,
      compound_id = row$compound_id,
      compound_name = row$compound_name,
      compound_name_clean = row$compound_name_clean,
      source_database = row$source_database,
      source_record_id = row$source_record_id,
      evidence_tier = row$evidence_tier,
      matched_rank = row$matched_rank,
      occurrence_status = row$occurrence_status,
      confidence = row$confidence,
      evidence_grade = grade$grade,
      evidence_grade_rank = grade$rank,
      evidence_grade_label = grade$label,
      evidence_grade_basis = grade$basis,
      source_backed = .bundle_yes_no(source_backed),
      structure_resolved = .bundle_yes_no(structure_resolved),
      comparable_for_analysis = .bundle_yes_no(comparable),
      review_required = .bundle_yes_no(review_required),
      curation_decision = decision,
      analysis_filter_key = paste(
        grade$grade,
        ifelse(structure_resolved, "structure_resolved",
               "structure_unresolved"),
        ifelse(comparable, "comparable", "not_comparable"),
        sep = "__"
      ),
      recommended_use = .evidence_grade_use_note(grade$grade,
                                                 structure_resolved,
                                                 comparable),
      stringsAsFactors = FALSE
    )
  })
  out = .bundle_bind(rows)
  out = out[, .evidence_grade_cols(), drop = FALSE]
  row.names(out) = NULL
  out
}

#' Summarize plant-pair Tanimoto by comparable chemistry
#'
#' @description
#' `plantComparableTanimotoSummary()` summarizes plant-compound pair Tanimoto
#' rows after requiring both compounds in a pair to share the same comparable
#' chemistry scope or group. This keeps like-with-like comparisons separate
#' from the unfiltered whole-chemistry plant-pair summary.
#'
#' @param plant_pair_tanimoto Optional unfiltered plant-pair summary.
#' @param plant_compound_pair_tanimoto Data frame or CSV/CSV.GZ path containing
#' plant-compound pair Tanimoto rows.
#' @param membership Enriched plant-compound membership table with
#' `comparison_scope`, `comparison_group`, and `comparable_for_matrix`.
#' @param include_unknown Logical. If `FALSE`, unknown/non-comparable chemistry
#' is excluded from filtered summaries.
#' @param thresholds Numeric Tanimoto thresholds to count.
#' @param max_pair_rows Maximum plant-compound pair rows to read from a file.
#'
#' @return List with `Overall`, `ScopeFiltered`, `GroupFiltered`, and
#' `Diagnostics`.
#'
#' @export
plantComparableTanimotoSummary = function(plant_pair_tanimoto = NULL,
                                          plant_compound_pair_tanimoto = NULL,
                                          membership,
                                          include_unknown = FALSE,
                                          thresholds = c(0.5, 0.7, 0.85,
                                                         0.95),
                                          max_pair_rows = Inf) {
  membership = .comparable_membership(membership, include_unknown)
  overall = if (is.data.frame(plant_pair_tanimoto)) {
    .bundle_enhance_pair_tanimoto(plant_pair_tanimoto)
  } else {
    .bundle_empty(.comparable_pair_summary_cols("overall"))
  }
  pairs = .read_optional_pair_table(plant_compound_pair_tanimoto,
                                    max_pair_rows = max_pair_rows)
  if (!is.data.frame(pairs) || nrow(pairs) < 1 || nrow(membership) < 1) {
    diag = data.frame(
      status = "not_computed",
      reason = "Comparable Tanimoto requires plant-compound pair rows and enriched comparable membership.",
      input_pair_rows = if (is.data.frame(pairs)) nrow(pairs) else 0L,
      comparable_membership_rows = nrow(membership),
      stringsAsFactors = FALSE
    )
    return(list(Overall = overall,
                ScopeFiltered = .bundle_empty(
                  .comparable_pair_summary_cols("scope", thresholds)
                ),
                GroupFiltered = .bundle_empty(
                  .comparable_pair_summary_cols("group", thresholds)
                ),
                Diagnostics = diag))
  }
  pairs = .normalize_plant_compound_pair_table(pairs)
  scope = .summarize_comparable_pairs(pairs, membership, "comparison_scope",
                                      "scope", thresholds)
  group = .summarize_comparable_pairs(pairs, membership, "comparison_group",
                                      "group", thresholds)
  diag = data.frame(
    status = "computed",
    reason = NA_character_,
    input_pair_rows = nrow(pairs),
    comparable_membership_rows = nrow(membership),
    scope_summary_rows = nrow(scope),
    group_summary_rows = nrow(group),
    stringsAsFactors = FALSE
  )
  list(Overall = overall, ScopeFiltered = scope, GroupFiltered = group,
       Diagnostics = diag)
}

#' Export model-ready plant chemistry feature tables
#'
#' @description
#' `exportPlantChemistryFeatureSet()` builds analysis-neutral feature tables
#' from enriched plant-compound membership rows. The tables are suitable for
#' joining to phylogenetic, ecological, water-quality, remediation, or other
#' project metadata.
#'
#' @param membership Enriched plant-compound membership table.
#' @param path Optional output directory. If supplied, CSV files and a manifest
#' are written.
#' @param overwrite Logical. If `TRUE`, replace an existing output directory.
#' @param modes Matrix modes to export. Supported values are `"count"`,
#' `"binary"`, `"fraction"`, and `"confidence"`. The default preserves the
#' original count-matrix output.
#' @param include_context Logical. If `TRUE`, include plant-part, tissue, and
#' method matrices when those fields are present.
#' @param include_evidence Logical. If `TRUE`, include evidence-grade matrices.
#' @param species_universe Optional character vector or data frame defining the
#' complete species set and row order for every feature matrix. Species present
#' in `membership` but absent from this input are appended. The default uses all
#' species in `membership`.
#' @param plant_metadata Optional plant metadata table used to populate
#' taxonomy fields for species with no membership rows.
#'
#' @return Named list of feature tables and manifest.
#'
#' @export
exportPlantChemistryFeatureSet = function(membership,
                                          path = NULL,
                                          overwrite = FALSE,
                                          modes = "count",
                                          include_context = TRUE,
                                          include_evidence = TRUE,
                                          species_universe = NULL,
                                          plant_metadata = NULL) {
  if (!is.data.frame(membership)) {
    stop("`membership` must be an enriched plant-compound membership table.",
         call. = FALSE)
  }
  membership = as.data.frame(membership, stringsAsFactors = FALSE)
  for (col in c("species", "compound_id", "comparison_group",
                "comparison_scope", "source_database",
                "comparable_for_matrix", "plant_part_group",
                "tissue_group", "method_group", "confidence")) {
    if (!col %in% names(membership)) membership[[col]] = NA_character_
  }
  species_universe = .feature_species_universe(membership, species_universe)
  if (!is.null(plant_metadata)) {
    plant_metadata = standardizePlantMetadata(plant_metadata)
  }
  modes = unique(tolower(.uaf_non_empty(modes)))
  modes = modes[modes %in% c("count", "binary", "fraction", "confidence")]
  if (length(modes) < 1) modes = "count"
  evidence = plantOccurrenceEvidenceGrade(membership)
  if (nrow(evidence) == nrow(membership)) {
    membership$evidence_grade = evidence$evidence_grade
  } else if (!"evidence_grade" %in% names(membership)) {
    membership$evidence_grade = NA_character_
  }
  comparable = membership[
    .bundle_truthy(membership$comparable_for_matrix) &
      !tolower(.bundle_squish(membership$comparison_group)) %in%
      c("", "unknown") &
      !tolower(.bundle_squish(membership$comparison_scope)) %in%
      c("", "unknown"),
    , drop = FALSE
  ]
  out = list()
  manifest_rows = list()
  add_matrix = function(name, table, field, mode, caveat) {
    out[[name]] <<- table
    manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
      Table = name,
      FileName = paste0(.feature_file_name(name), ".csv"),
      RowCount = nrow(table),
      ColumnCount = ncol(table),
      FeatureField = field,
      Mode = mode,
      SpeciesUniverseCount = length(species_universe),
      EvidenceFilter = if (field %in% c("comparison_group",
                                        "comparison_scope")) {
        "comparable_for_matrix == Yes and non-unknown scope/group"
      } else {
        "all supplied membership rows"
      },
      Caveat = caveat,
      stringsAsFactors = FALSE
    )
  }
  mode_label = function(x) paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  for (mode in modes) {
    suffix = mode_label(mode)
    add_matrix(paste0("ComparisonGroup", suffix, "Matrix"),
               .feature_matrix(comparable, "comparison_group", mode,
                               species_universe),
               "comparison_group", mode,
               "Comparable chemistry group features are source-backed or explicitly classified; unknown chemistry is excluded.")
    add_matrix(paste0("ComparisonScope", suffix, "Matrix"),
               .feature_matrix(comparable, "comparison_scope", mode,
                               species_universe),
               "comparison_scope", mode,
               "Comparable chemistry scope features separate unlike chemistry before downstream modeling.")
    add_matrix(paste0("SourceCoverage", suffix, "Matrix"),
               .feature_matrix(membership, "source_database", mode,
                               species_universe),
               "source_database", mode,
               "Source coverage reflects public/source records, not biological completeness.")
    if (isTRUE(include_evidence)) {
      add_matrix(paste0("EvidenceGrade", suffix, "Matrix"),
                 .feature_matrix(membership, "evidence_grade", mode,
                                 species_universe),
                 "evidence_grade", mode,
                 "Evidence grades are conservative analysis tiers, not experimental confirmation.")
    }
    if (isTRUE(include_context)) {
      add_matrix(paste0("PlantPart", suffix, "Matrix"),
                 .feature_matrix(membership, "plant_part_group", mode,
                                 species_universe),
                 "plant_part_group", mode,
                 "Plant-part context is included only when reported or extracted from sources.")
      add_matrix(paste0("Tissue", suffix, "Matrix"),
                 .feature_matrix(membership, "tissue_group", mode,
                                 species_universe),
                 "tissue_group", mode,
                 "Tissue context is included only when reported or extracted from sources.")
      add_matrix(paste0("Method", suffix, "Matrix"),
                 .feature_matrix(membership, "method_group", mode,
                                 species_universe),
                 "method_group", mode,
                 "Method context reflects source metadata and is often incomplete.")
    }
  }
  group_count = out$ComparisonGroupCountMatrix
  scope_count = out$ComparisonScopeCountMatrix
  source_count = out$SourceCoverageCountMatrix
  metadata = .feature_species_metadata(membership, species_universe,
                                       plant_metadata)
  out$SourceCoverageMatrix = source_count
  out$SpeciesMetadata = metadata
  manifest_rows[[length(manifest_rows) + 1L]] = data.frame(
    Table = "SpeciesMetadata",
    FileName = "species_metadata.csv",
    RowCount = nrow(metadata),
    ColumnCount = ncol(metadata),
    FeatureField = "species_quality_metadata",
    Mode = "metadata",
    SpeciesUniverseCount = length(species_universe),
    EvidenceFilter = "all supplied membership rows",
    Caveat = "Quality metadata summarize evidence coverage and should not be interpreted as biological completeness.",
    stringsAsFactors = FALSE
  )
  manifest = .bundle_bind(manifest_rows)
  out$Manifest = manifest
  if (!is.null(path)) {
    if ((dir.exists(path) || file.exists(path)) && !isTRUE(overwrite)) {
      stop("Feature-set output already exists. Use `overwrite = TRUE`: ",
           path, call. = FALSE)
    }
    if (dir.exists(path) || file.exists(path)) {
      unlink(path, recursive = TRUE, force = TRUE)
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    for (i in seq_len(nrow(manifest))) {
      .categorate_write_csv_file(out[[manifest$Table[[i]]]],
                                 file.path(path, manifest$FileName[[i]]))
    }
    .categorate_write_csv_file(manifest,
                               file.path(path, "feature_manifest.csv"))
  }
  out
}

#' Run a reusable plant chemistry project workflow
#'
#' @description
#' `runPlantChemistryProject()` coordinates plant list/metadata intake,
#' curated or discovered plant-compound rows, optional cached categorate batch
#' exports, bundle finalization, validation, evidence grading, and feature-set
#' outputs. Live discovery is opt-in through `run_discovery = TRUE`.
#'
#' @param plant_list Character vector, data frame, or CSV path of plant names.
#' @param output_dir Project output directory.
#' @param metadata Optional plant metadata table or path.
#' @param plant_compounds Optional curated plant-compound table or path.
#' @param plant_result Optional result object or RDS path from plant workflows.
#' @param categorate_batches Optional categorate batch directory or files.
#' @param species_pair_tanimoto Optional plant-pair Tanimoto summary table/path.
#' @param resolved_compounds Optional resolved compounds table/path.
#' @param pubchem_fingerprints Optional PubChem fingerprints table/path.
#' @param plant_compound_pair_tanimoto Optional plant-compound pair table/path
#' for comparable Tanimoto summaries.
#' @param file_references Optional file references recorded in the bundle.
#' @param project_id Optional project label.
#' @param run_discovery Logical. If `TRUE`, call `resolvePlantPhytochemistry()`.
#' @param overwrite Logical. If `TRUE`, replace existing project output.
#' @param ... Additional arguments passed to `resolvePlantPhytochemistry()` when
#' `run_discovery = TRUE`.
#'
#' @return Project manifest list.
#'
#' @export
runPlantChemistryProject = function(plant_list,
                                    output_dir,
                                    metadata = NULL,
                                    plant_compounds = NULL,
                                    plant_result = NULL,
                                    categorate_batches = NULL,
                                    species_pair_tanimoto = NULL,
                                    resolved_compounds = NULL,
                                    pubchem_fingerprints = NULL,
                                    plant_compound_pair_tanimoto = NULL,
                                    file_references = NULL,
                                    project_id = NULL,
                                    run_discovery = FALSE,
                                    overwrite = FALSE,
                                    ...) {
  if (missing(output_dir) || length(.uaf_non_empty(output_dir)) != 1) {
    stop("`output_dir` is required.", call. = FALSE)
  }
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if ((dir.exists(output_dir) || file.exists(output_dir)) && !isTRUE(overwrite)) {
    stop("Project output exists. Use `overwrite = TRUE`: ", output_dir,
         call. = FALSE)
  }
  if (dir.exists(output_dir) || file.exists(output_dir)) {
    unlink(output_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metadata_std = if (!is.null(metadata)) standardizePlantMetadata(metadata) else
    NULL
  plants = .bundle_standardize_plant_list(plant_list, metadata_std)
  if (nrow(plants) < 1) stop("No plant names were supplied.", call. = FALSE)

  result = .project_read_result(plant_result)
  if (is.null(result) && isTRUE(run_discovery)) {
    result = resolvePlantPhytochemistry(plants = plants$species, ...)
    saveRDS(result, file.path(output_dir, "plant_phytochemistry_result.rds"))
  }
  membership = .project_membership_from_inputs(result, plant_compounds)
  if (!is.data.frame(membership) || nrow(membership) < 1) {
    stop("Supply `plant_compounds`, `plant_result`, or set ",
         "`run_discovery = TRUE`.", call. = FALSE)
  }
  membership_file = file.path(output_dir,
                              "project_plant_compound_membership.csv")
  .categorate_write_csv_file(membership, membership_file)
  if (!is.null(metadata_std)) {
    .categorate_write_csv_file(metadata_std,
                               file.path(output_dir,
                                         "project_plant_metadata.csv"))
  }

  bundle_dir = file.path(output_dir, "plant_chemistry_analysis_bundle")
  if (!is.null(categorate_batches)) {
    bundle_manifest = exportPlantChemistryAnalysisBundle(
      categorate_batches = categorate_batches,
      path = bundle_dir,
      plant_membership = membership_file,
      species_pair_tanimoto = species_pair_tanimoto,
      resolved_compounds = resolved_compounds,
      pubchem_fingerprints = pubchem_fingerprints,
      plant_compound_pair_tanimoto = plant_compound_pair_tanimoto,
      file_references = file_references,
      plant_list = plants,
      metadata = metadata_std,
      project_id = project_id,
      overwrite = TRUE,
      finalize = TRUE,
      include_feature_exports = TRUE,
      include_comparable_tanimoto = !is.null(plant_compound_pair_tanimoto),
      validate_export = TRUE
    )
  } else {
    .project_write_minimal_bundle(bundle_dir, membership, plants, metadata_std,
                                  project_id = project_id)
    bundle_manifest = finalizePlantChemistryAnalysisBundle(
      bundle_dir,
      plant_list = plants,
      metadata = metadata_std,
      plant_compound_pair_tanimoto = plant_compound_pair_tanimoto,
      project_id = project_id,
      overwrite = TRUE,
      include_comparable_tanimoto = !is.null(plant_compound_pair_tanimoto),
      validate_export = TRUE
    )
  }
  validation = validatePlantChemistryAnalysisBundle(bundle_dir)
  manifest = list(
    Project = data.frame(
      project_id = .uaf_first_non_empty_text(project_id, "plant_chemistry_project"),
      output_dir = output_dir,
      bundle_dir = bundle_dir,
      run_discovery = isTRUE(run_discovery),
      bundle_validation_status = validation$Summary$ExportReadyStatus[[1]],
      created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      stringsAsFactors = FALSE
    ),
    Inputs = data.frame(
      input = c("plant_list", "metadata", "plant_compounds",
                "plant_result", "categorate_batches", "species_pair_tanimoto",
                "resolved_compounds", "pubchem_fingerprints",
                "plant_compound_pair_tanimoto"),
      supplied = c(TRUE, !is.null(metadata), !is.null(plant_compounds),
                   !is.null(plant_result), !is.null(categorate_batches),
                   !is.null(species_pair_tanimoto), !is.null(resolved_compounds),
                   !is.null(pubchem_fingerprints),
                   !is.null(plant_compound_pair_tanimoto)),
      stringsAsFactors = FALSE
    ),
    BundleManifest = bundle_manifest,
    ValidationSummary = validation$Summary
  )
  .categorate_write_csv_file(manifest$Project,
                             file.path(output_dir, "project_run_manifest.csv"))
  jsonlite::write_json(manifest,
                       file.path(output_dir, "project_run_manifest.json"),
                       pretty = TRUE, auto_unbox = TRUE, null = "null")
  class(manifest) = c("uaf_plant_chemistry_project", "list")
  manifest
}

.metadata_detect_col = function(x, requested = NULL, candidates = character()) {
  if (!is.data.frame(x) || ncol(x) < 1) return(NULL)
  nm = names(x)
  if (!is.null(requested) && length(requested) == 1 &&
      !is.na(requested) && nzchar(requested)) {
    if (requested %in% nm) return(requested)
    requested_norm = .plant_normalize_column_names(requested)
    if (requested_norm %in% nm) return(requested_norm)
  }
  candidates = .plant_normalize_column_names(candidates)
  hit = candidates[candidates %in% nm]
  if (length(hit) > 0) hit[[1]] else NULL
}

.evidence_grade_cols = function() {
  c("row_id", "species", "compound_id", "compound_name",
    "compound_name_clean", "source_database", "source_record_id",
    "evidence_tier", "matched_rank", "occurrence_status", "confidence",
    "evidence_grade", "evidence_grade_rank", "evidence_grade_label",
    "evidence_grade_basis", "source_backed", "structure_resolved",
    "comparable_for_analysis", "review_required", "curation_decision",
    "analysis_filter_key", "recommended_use")
}

.evidence_review_decisions = function(occurrences, review_table = NULL) {
  decision = rep("not_reviewed", nrow(occurrences))
  if (is.null(review_table)) {
    return(data.frame(review_decision = decision, stringsAsFactors = FALSE))
  }
  if (is.character(review_table) && length(review_table) == 1 &&
      file.exists(review_table)) {
    review_table = utils::read.csv(review_table, stringsAsFactors = FALSE,
                                   check.names = FALSE)
  }
  if (!is.data.frame(review_table) || nrow(review_table) < 1) {
    return(data.frame(review_decision = decision, stringsAsFactors = FALSE))
  }
  review = as.data.frame(review_table, stringsAsFactors = FALSE)
  names(review) = .plant_normalize_column_names(names(review))
  decision_col = .metadata_detect_col(
    review, NULL,
    c("curation_decision", "review_decision", "decision", "status")
  )
  if (is.null(decision_col)) {
    return(data.frame(review_decision = decision, stringsAsFactors = FALSE))
  }
  if ("row_id" %in% names(review)) {
    idx = suppressWarnings(as.integer(review$row_id))
    ok = !is.na(idx) & idx >= 1 & idx <= length(decision)
    decision[idx[ok]] = .bundle_first_non_empty(review[[decision_col]][ok],
                                                decision[idx[ok]])
  }
  occ_key = .bundle_occurrence_key(occurrences)
  if ("occurrence_key" %in% names(review)) {
    idx = match(occ_key, review$occurrence_key)
  } else if (all(c("species", "compound_name_clean", "source_database",
                   "source_record_id") %in% names(review))) {
    idx = match(occ_key, .bundle_occurrence_key(review))
  } else if (all(c("species", "compound_name") %in% names(review))) {
    left = paste(occurrences$species, occurrences$compound_name, sep = "\r")
    right = paste(review$species, review$compound_name, sep = "\r")
    idx = match(left, right)
  } else {
    idx = rep(NA_integer_, length(decision))
  }
  hit = !is.na(idx)
  decision[hit] = .bundle_first_non_empty(review[[decision_col]][idx[hit]],
                                          decision[hit])
  data.frame(review_decision = .bundle_squish(decision),
             stringsAsFactors = FALSE)
}

.evidence_grade_for_row = function(row, decision = "not_reviewed") {
  decision_clean = tolower(.bundle_squish(decision))
  if (decision_clean %in% c("exclude", "excluded", "remove",
                            "do_not_use", "reject")) {
    return(list(grade = "excluded_by_review", rank = 99L,
                label = "Excluded by review",
                basis = "A review table marked this occurrence as excluded."))
  }

  tier = tolower(.bundle_squish(row$evidence_tier))
  source = tolower(.bundle_squish(row$source_database))
  rank = tolower(.bundle_squish(row$matched_rank))
  status = tolower(.bundle_squish(row$occurrence_status))
  has_source_record = .bundle_known(row$source_record_id) |
    .bundle_known(row$evidence_url) | .bundle_known(row$pmid) |
    .bundle_known(row$doi)
  has_compound = .bundle_known(row$compound_name) |
    .bundle_known(row$compound_name_clean) | .bundle_known(row$compound_id)

  if (decision_clean %in% c("promote_direct_species_database",
                            "direct_species_database",
                            "confirmed_direct_database")) {
    return(list(grade = "direct_species_database_record", rank = 1L,
                label = "Direct species database record",
                basis = "A review table promoted this row to direct species database evidence."))
  }
  if (decision_clean %in% c("promote_direct_species_literature",
                            "direct_species_literature",
                            "confirmed_literature")) {
    return(list(grade = "direct_species_literature_supported_record",
                rank = 2L,
                label = "Direct species literature-supported record",
                basis = "A review table promoted this row to direct species literature evidence."))
  }
  if (!has_compound || grepl("unresolved", tier) ||
      grepl("unresolved", status)) {
    return(list(grade = "unresolved_or_review_required", rank = 5L,
                label = "Unresolved or review-required",
                basis = "The occurrence lacks a resolved compound or is explicitly unresolved."))
  }
  if (grepl("direct_species_database", tier) ||
      (identical(rank, "species") &&
       !grepl("pubmed|pubtator|literature|candidate", source) &&
       has_source_record)) {
    return(list(grade = "direct_species_database_record", rank = 1L,
                label = "Direct species database record",
                basis = "The row is labeled as species-level database evidence or has a species-level source record."))
  }
  if (grepl("direct_species_literature", tier) ||
      (identical(rank, "species") &&
       (grepl("pubmed|literature", source) | .bundle_known(row$pmid) |
          .bundle_known(row$doi)) &&
       !grepl("candidate|pubtator", tier))) {
    return(list(grade = "direct_species_literature_supported_record",
                rank = 2L,
                label = "Direct species literature-supported record",
                basis = "The row links a species-level occurrence to literature evidence."))
  }
  if (grepl("genus|family|fallback", tier) || rank %in% c("genus", "family")) {
    return(list(grade = "source_backed_genus_family_fallback", rank = 3L,
                label = "Source-backed genus/family fallback",
                basis = "The record is source-backed but not direct species-level evidence."))
  }
  if (grepl("pubtator|candidate|co.?mention", tier) ||
      grepl("pubtator", source) ||
      (grepl("pubmed", source) && !grepl("direct", tier))) {
    return(list(grade = "pubtator_pubmed_candidate_only", rank = 4L,
                label = "PubTator/PubMed candidate only",
                basis = "The row is literature-mining or co-mention evidence and should be reviewed before use as occurrence evidence."))
  }
  list(grade = "unresolved_or_review_required", rank = 5L,
       label = "Unresolved or review-required",
       basis = "The row does not meet a stronger direct or fallback evidence rule.")
}

.evidence_grade_use_note = function(grade, structure_resolved, comparable) {
  if (identical(grade, "excluded_by_review")) {
    return("Do not use in downstream analyses unless the review decision is changed.")
  }
  if (identical(grade, "pubtator_pubmed_candidate_only")) {
    return("Use only for literature triage until a source-backed occurrence is curated.")
  }
  if (identical(grade, "source_backed_genus_family_fallback")) {
    return("Use as broader taxonomic context; keep separate from direct species evidence.")
  }
  if (identical(grade, "unresolved_or_review_required")) {
    return("Review before analysis; do not treat as a resolved plant-compound occurrence.")
  }
  if (!isTRUE(structure_resolved)) {
    return("Useful for occurrence summaries but not structure-based Tanimoto analyses.")
  }
  if (!isTRUE(comparable)) {
    return("Useful for source-backed summaries; exclude from comparable chemistry matrices.")
  }
  "Use in source-backed comparable chemistry analyses with public-data caveats."
}

.bundle_evidence_grade_summary = function(evidence_grades) {
  cols = c("evidence_grade", "evidence_grade_rank", "evidence_grade_label",
           "occurrence_count", "species_count", "compound_count",
           "source_backed_count", "structure_resolved_count",
           "comparable_count", "review_required_count", "recommended_use")
  if (!is.data.frame(evidence_grades) || nrow(evidence_grades) < 1) {
    return(.bundle_empty(cols))
  }
  rows = lapply(split(evidence_grades, evidence_grades$evidence_grade),
                function(x) {
    data.frame(
      evidence_grade = x$evidence_grade[[1]],
      evidence_grade_rank =
        suppressWarnings(as.integer(x$evidence_grade_rank[[1]])),
      evidence_grade_label = x$evidence_grade_label[[1]],
      occurrence_count = nrow(x),
      species_count = length(unique(.uaf_non_empty(x$species))),
      compound_count = length(unique(.uaf_non_empty(x$compound_id))),
      source_backed_count = sum(.bundle_truthy(x$source_backed), na.rm = TRUE),
      structure_resolved_count = sum(.bundle_truthy(x$structure_resolved),
                                     na.rm = TRUE),
      comparable_count = sum(.bundle_truthy(x$comparable_for_analysis),
                             na.rm = TRUE),
      review_required_count = sum(.bundle_truthy(x$review_required),
                                  na.rm = TRUE),
      recommended_use = x$recommended_use[[1]],
      stringsAsFactors = FALSE
    )
  })
  out = .bundle_bind(rows)
  out = out[order(out$evidence_grade_rank, out$evidence_grade), , drop = FALSE]
  row.names(out) = NULL
  out[, cols, drop = FALSE]
}

.bundle_review_required_occurrences = function(enriched, evidence_grades) {
  cols = c("row_id", "species", "compound_id", "compound_name",
           "source_database", "source_record_id", "evidence_tier",
           "matched_rank", "evidence_grade", "review_category",
           "review_reason", "recommended_action",
           "evidence_review_required", "identity_review_required",
           "structure_review_required", "context_review_required",
           "comparability_review_required", "citation_review_required",
           "recommended_use", "evidence_url", "pmid", "doi")
  if (!is.data.frame(evidence_grades) || nrow(evidence_grades) < 1) {
    return(.bundle_empty(cols))
  }
  enriched = as.data.frame(enriched, stringsAsFactors = FALSE)
  for (col in c("evidence_url", "pmid", "doi", "identity_review_required",
                "identity_ambiguity_flag", "identity_review_reason",
                "identity_recommended_action", "biological_context_known",
                "context_known_record", "comparison_scope")) {
    if (!col %in% names(enriched)) enriched[[col]] = NA_character_
  }
  all_idx = suppressWarnings(as.integer(evidence_grades$row_id))
  identity_flag = .bundle_truthy(.bundle_index(
    enriched, all_idx, "identity_review_required"
  )) | .bundle_truthy(.bundle_index(
    enriched, all_idx, "identity_ambiguity_flag"
  ))
  structure_flag = !.bundle_truthy(evidence_grades$structure_resolved)
  biological_context = .bundle_first_non_empty(
    .bundle_index(enriched, all_idx, "biological_context_known"),
    .bundle_index(enriched, all_idx, "context_known_record")
  )
  context_flag = !.bundle_truthy(biological_context)
  scope = tolower(.bundle_squish(.bundle_index(
    enriched, all_idx, "comparison_scope"
  )))
  comparability_flag = is.na(scope) | scope == "" | scope %in%
    c("unknown", "broad_or_uncertain", "xenobiotic_or_contaminant")
  evidence_flag = .bundle_truthy(evidence_grades$review_required) |
    evidence_grades$evidence_grade %in%
    c("unresolved_or_review_required", "pubtator_pubmed_candidate_only",
      "source_backed_genus_family_fallback", "excluded_by_review")
  has_citation = .bundle_known(.bundle_index(enriched, all_idx,
                                             "source_record_id")) |
    .bundle_known(.bundle_index(enriched, all_idx, "evidence_url")) |
    .bundle_known(.bundle_index(enriched, all_idx, "pmid")) |
    .bundle_known(.bundle_index(enriched, all_idx, "doi"))
  citation_flag = evidence_grades$source_backed == "Yes" & !has_citation
  keep = evidence_flag | identity_flag | structure_flag | context_flag |
    comparability_flag | citation_flag
  review = evidence_grades[keep, , drop = FALSE]
  if (nrow(review) < 1) return(.bundle_empty(cols))
  identity_flag = identity_flag[keep]
  structure_flag = structure_flag[keep]
  context_flag = context_flag[keep]
  comparability_flag = comparability_flag[keep]
  evidence_flag = evidence_flag[keep]
  citation_flag = citation_flag[keep]
  idx = all_idx[keep]
  details = lapply(seq_len(nrow(review)), function(i) {
    categories = reasons = actions = character()
    if (evidence_flag[[i]]) {
      categories = c(categories, "evidence")
      reasons = c(reasons, if (review$evidence_grade[[i]] %in%
                                c("direct_species_database_record",
                                  "direct_species_literature_supported_record")) {
        "A source curation flag recommends checking the species-compound record."
      } else {
        review$evidence_grade_basis[[i]]
      })
      actions = c(actions, "Verify the source record and retain the original evidence grade unless stronger evidence is documented.")
    }
    if (identity_flag[[i]]) {
      categories = c(categories, "identity")
      reasons = c(reasons, .uaf_first_non_empty_text(
        .bundle_index(enriched, idx[[i]], "identity_review_reason"),
        "The compound identity is ambiguous or explicitly marked for review."
      ))
      actions = c(actions, .uaf_first_non_empty_text(
        .bundle_index(enriched, idx[[i]], "identity_recommended_action"),
        "Resolve, replace, exclude, or caveat the identity in a reproducible review table."
      ))
    }
    if (structure_flag[[i]]) {
      categories = c(categories, "structure")
      reasons = c(reasons, "No usable SMILES, InChIKey, or PubChem CID is available for structure-based analysis.")
      actions = c(actions, "Exclude from Tanimoto analysis unless a source-backed structure is resolved.")
    }
    if (context_flag[[i]]) {
      categories = c(categories, "biological_context")
      reasons = c(reasons, "No source-backed plant-part or tissue context is available.")
      actions = c(actions, "Exclude from plant-part or tissue-specific comparisons unless context is curated from a source.")
    }
    if (comparability_flag[[i]]) {
      categories = c(categories, "comparability")
      reasons = c(reasons, "The chemistry comparison scope is unknown, broad/uncertain, or non-comparable by default.")
      actions = c(actions, "Exclude from comparable chemistry matrices unless a source-backed classification override is recorded.")
    }
    if (citation_flag[[i]]) {
      categories = c(categories, "citation")
      reasons = c(reasons, "The source-backed grade lacks a record identifier, URL, PMID, or DOI in this bundle.")
      actions = c(actions, "Add a traceable source reference or downgrade the evidence grade.")
    }
    data.frame(
      review_category = paste(unique(categories), collapse = "; "),
      review_reason = paste(unique(.uaf_non_empty(reasons)), collapse = "; "),
      recommended_action = paste(unique(.uaf_non_empty(actions)),
                                 collapse = "; "),
      stringsAsFactors = FALSE
    )
  })
  details = .bundle_bind(details)
  out = data.frame(
    row_id = review$row_id,
    species = review$species,
    compound_id = review$compound_id,
    compound_name = review$compound_name,
    source_database = review$source_database,
    source_record_id = review$source_record_id,
    evidence_tier = review$evidence_tier,
    matched_rank = review$matched_rank,
    evidence_grade = review$evidence_grade,
    review_category = details$review_category,
    review_reason = details$review_reason,
    recommended_action = details$recommended_action,
    evidence_review_required = .bundle_yes_no(evidence_flag),
    identity_review_required = .bundle_yes_no(identity_flag),
    structure_review_required = .bundle_yes_no(structure_flag),
    context_review_required = .bundle_yes_no(context_flag),
    comparability_review_required = .bundle_yes_no(comparability_flag),
    citation_review_required = .bundle_yes_no(citation_flag),
    recommended_use = ifelse(
      identity_flag | structure_flag | comparability_flag |
        review$evidence_grade %in%
        c("unresolved_or_review_required", "pubtator_pubmed_candidate_only",
          "excluded_by_review"),
      "Retain for audit; exclude from the affected analysis until review is completed.",
      ifelse(context_flag,
             "May support non-contextual occurrence summaries; exclude from plant-part or tissue-specific analysis until context is resolved.",
             review$recommended_use)
    ),
    evidence_url = .bundle_index(enriched, idx, "evidence_url"),
    pmid = .bundle_index(enriched, idx, "pmid"),
    doi = .bundle_index(enriched, idx, "doi"),
    stringsAsFactors = FALSE
  )
  out[, cols, drop = FALSE]
}

.comparable_membership = function(membership, include_unknown = FALSE) {
  if (!is.data.frame(membership) || nrow(membership) < 1) {
    return(data.frame())
  }
  x = as.data.frame(membership, stringsAsFactors = FALSE)
  for (col in c("species", "compound_id", "comparison_scope",
                "comparison_group", "comparable_for_matrix",
                "evidence_tier", "source_database")) {
    if (!col %in% names(x)) x[[col]] = NA_character_
  }
  x$comparison_scope = .bundle_first_non_empty(x$comparison_scope, "unknown")
  x$comparison_group = .bundle_first_non_empty(x$comparison_group, "unknown")
  ok = .bundle_truthy(x$comparable_for_matrix)
  if (!isTRUE(include_unknown)) {
    ok = ok &
      !tolower(.bundle_squish(x$comparison_scope)) %in% c("", "unknown") &
      !tolower(.bundle_squish(x$comparison_group)) %in% c("", "unknown")
  }
  x = x[ok, , drop = FALSE]
  x = x[.bundle_known(x$species) & .bundle_known(x$compound_id), ,
        drop = FALSE]
  key = paste(x$species, x$compound_id, x$comparison_scope,
              x$comparison_group, sep = "\r")
  x = x[!duplicated(key), , drop = FALSE]
  row.names(x) = NULL
  x
}

.read_optional_pair_table = function(x, max_pair_rows = Inf) {
  if (is.null(x)) return(data.frame())
  if (is.data.frame(x)) return(as.data.frame(x, stringsAsFactors = FALSE))
  path = .uaf_non_empty(x)
  if (length(path) != 1 || !file.exists(path)) {
    stop("Plant-compound pair table path does not exist.", call. = FALSE)
  }
  nrows = if (is.finite(max_pair_rows)) as.integer(max_pair_rows) else -1L
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                  nrows = nrows)
}

.normalize_plant_compound_pair_table = function(pairs) {
  pairs = as.data.frame(pairs, stringsAsFactors = FALSE)
  names(pairs) = .plant_normalize_column_names(names(pairs))
  rename = c(group_a = "species_a", group_b = "species_b",
             query_a = "species_a", query_b = "species_b")
  for (from in names(rename)) {
    to = rename[[from]]
    if (!to %in% names(pairs) && from %in% names(pairs)) {
      names(pairs)[names(pairs) == from] = to
    }
  }
  for (col in c("species_a", "species_b", "compound_id_a",
                "compound_id_b", "compound_name_a", "compound_name_b",
                "tanimoto")) {
    if (!col %in% names(pairs)) pairs[[col]] = NA_character_
  }
  pairs$tanimoto = suppressWarnings(as.numeric(pairs$tanimoto))
  pairs = pairs[.bundle_known(pairs$species_a) &
                  .bundle_known(pairs$species_b) &
                  .bundle_known(pairs$compound_id_a) &
                  .bundle_known(pairs$compound_id_b) &
                  !is.na(pairs$tanimoto), , drop = FALSE]
  row.names(pairs) = NULL
  pairs
}

.summarize_comparable_pairs = function(pairs, membership, field, type,
                                       thresholds) {
  cols = .comparable_pair_summary_cols(type, thresholds)
  if (nrow(pairs) < 1 || nrow(membership) < 1 || !field %in% names(membership)) {
    return(.bundle_empty(cols))
  }
  left_key = paste(pairs$species_a, pairs$compound_id_a, sep = "\r")
  right_key = paste(pairs$species_b, pairs$compound_id_b, sep = "\r")
  membership_key = paste(membership$species, membership$compound_id, sep = "\r")
  ia = match(left_key, membership_key)
  ib = match(right_key, membership_key)
  value_a = .bundle_index(membership, ia, field)
  value_b = .bundle_index(membership, ib, field)
  ok = !is.na(ia) & !is.na(ib) & .bundle_known(value_a) &
    .bundle_known(value_b) & value_a == value_b
  if (!any(ok)) return(.bundle_empty(cols))
  x = pairs[ok, , drop = FALSE]
  x$comparison_value = value_a[ok]
  x$species_pair_key = paste(pmin(x$species_a, x$species_b),
                             pmax(x$species_a, x$species_b), sep = "\r")
  x$summary_key = paste(x$species_pair_key, x$comparison_value, sep = "\r")
  rows = lapply(split(x, x$summary_key), function(g) {
    tan = suppressWarnings(as.numeric(g$tanimoto))
    tan = tan[!is.na(tan)]
    if (length(tan) < 1) return(NULL)
    species_a = pmin(g$species_a[[1]], g$species_b[[1]])
    species_b = pmax(g$species_a[[1]], g$species_b[[1]])
    count_cols = stats::setNames(
      as.list(vapply(thresholds, function(th) sum(tan >= th, na.rm = TRUE),
                     integer(1))),
      paste0("compound_pair_count_ge_", gsub("[.]", "_", thresholds))
    )
    support = .comparable_support(length(tan))
    base = data.frame(
      species_a = species_a,
      species_b = species_b,
      unordered_pair_key = paste(species_a, species_b, sep = " || "),
      comparison_type = type,
      comparison_value = g$comparison_value[[1]],
      compound_pair_count = length(tan),
      mean_tanimoto = round(mean(tan), 6),
      median_tanimoto = round(stats::median(tan), 6),
      p95_tanimoto = round(as.numeric(stats::quantile(tan, 0.95,
                                                       names = FALSE,
                                                       na.rm = TRUE)), 6),
      max_tanimoto = round(max(tan), 6),
      support_tier = support$tier,
      low_support_caution = .bundle_yes_no(support$low),
      support_note = support$note,
      comparison_filter =
        paste0(field, " == '", g$comparison_value[[1]], "'"),
      caution =
        "Filtered summary includes only compound pairs with matching comparable chemistry classification.",
      stringsAsFactors = FALSE
    )
    if (identical(type, "scope")) {
      base$comparison_scope = g$comparison_value[[1]]
    } else {
      base$comparison_group = g$comparison_value[[1]]
    }
    cbind(base, as.data.frame(count_cols, stringsAsFactors = FALSE))
  })
  out = .bundle_bind(rows)
  for (col in setdiff(cols, names(out))) out[[col]] = NA
  out = out[, cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.comparable_support = function(n) {
  tier = if (is.na(n)) "unknown" else if (n >= 10000) {
    "high"
  } else if (n >= 1000) {
    "moderate"
  } else if (n >= 100) {
    "low"
  } else {
    "very_low"
  }
  list(tier = tier,
       low = !is.na(n) && n < 100,
       note = if (is.na(n)) {
         "Comparable compound-pair support was not reported."
       } else if (n < 100) {
         "Interpret with caution: fewer than 100 comparable compound-pair comparisons support this filtered summary."
       } else {
         "Filtered pair summary has at least 100 comparable compound-pair comparisons."
       })
}

.comparable_pair_summary_cols = function(type = c("scope", "group",
                                                  "overall"),
                                          thresholds = c(0.5, 0.7, 0.85,
                                                         0.95)) {
  type = match.arg(type)
  value_col = if (identical(type, "scope")) {
    "comparison_scope"
  } else if (identical(type, "group")) {
    "comparison_group"
  } else {
    character()
  }
  c("species_a", "species_b", "unordered_pair_key", "comparison_type",
    value_col, "comparison_value", "compound_pair_count", "mean_tanimoto",
    "median_tanimoto", "p95_tanimoto", "max_tanimoto",
    paste0("compound_pair_count_ge_", gsub("[.]", "_", thresholds)),
    "support_tier", "low_support_caution", "support_note",
    "comparison_filter", "caution")
}

.feature_matrix = function(x, field, mode = c("count", "binary",
                                              "fraction", "confidence"),
                            species_universe = NULL) {
  mode = match.arg(mode)
  cols = c("species_id", "species")
  if (!is.data.frame(x)) x = data.frame()
  species = .feature_species_universe(x, species_universe)
  if (length(species) < 1) return(.bundle_empty(cols))
  if (nrow(x) < 1 || !all(c("species", field) %in% names(x))) {
    return(data.frame(species_id = .feature_species_id(species),
                      species = species, stringsAsFactors = FALSE))
  }
  x = as.data.frame(x, stringsAsFactors = FALSE)
  x = x[.bundle_known(x$species) & .bundle_known(x[[field]]), , drop = FALSE]
  if (nrow(x) < 1) {
    return(data.frame(species_id = .feature_species_id(species),
                      species = species, stringsAsFactors = FALSE))
  }
  x$value = .feature_key(field, x[[field]])
  x$confidence_weight = if ("confidence" %in% names(x)) {
    .plant_confidence_score(x$confidence)
  } else {
    rep(1, nrow(x))
  }
  values = sort(unique(x$value))
  mat = matrix(0, nrow = length(species), ncol = length(values),
               dimnames = list(species, values))
  compound_total = stats::setNames(rep(0L, length(species)), species)
  for (sp in species) {
    hit = x[x$species == sp, , drop = FALSE]
    compound_total[[sp]] = length(unique(.uaf_non_empty(hit$compound_id)))
    for (value in values) {
      value_hit = hit[hit$value == value, , drop = FALSE]
      if (identical(mode, "confidence")) {
        compound_groups = split(value_hit$confidence_weight,
                                value_hit$compound_id)
        mat[sp, value] = round(sum(vapply(compound_groups, max, numeric(1)),
                                   na.rm = TRUE), 6)
      } else {
        mat[sp, value] = length(unique(.uaf_non_empty(value_hit$compound_id)))
      }
    }
  }
  if (identical(mode, "binary")) mat = ifelse(mat > 0, 1L, 0L)
  if (identical(mode, "fraction")) {
    denom = pmax(compound_total[rownames(mat)], 1L)
    mat = round(mat / denom, 6)
  }
  out = data.frame(species_id = .feature_species_id(species),
                   species = species,
                   as.data.frame(mat, check.names = FALSE),
                   stringsAsFactors = FALSE)
  row.names(out) = NULL
  out
}

.feature_file_name = function(name) {
  x = gsub("([a-z0-9])([A-Z])", "\\1_\\2", name)
  x = gsub("[^A-Za-z0-9]+", "_", x)
  tolower(gsub("^_+|_+$", "", x))
}

.bundle_finalize_feature_manifest = function(manifest) {
  if (!is.data.frame(manifest) || nrow(manifest) < 1) return(manifest)
  files = c(
    ComparisonGroupCountMatrix =
      "18_FeatureComparisonGroupCountMatrix.csv",
    ComparisonScopeCountMatrix =
      "19_FeatureComparisonScopeCountMatrix.csv",
    SourceCoverageCountMatrix = "20_FeatureSourceCoverageMatrix.csv",
    SpeciesMetadata = "21_FeatureSpeciesMetadata.csv",
    EvidenceGradeCountMatrix = "24_FeatureEvidenceGradeCountMatrix.csv",
    PlantPartCountMatrix = "25_FeaturePlantPartCountMatrix.csv",
    TissueCountMatrix = "26_FeatureTissueCountMatrix.csv",
    MethodCountMatrix = "27_FeatureMethodCountMatrix.csv"
  )
  idx = match(manifest$Table, names(files))
  hit = !is.na(idx)
  manifest$FileName[hit] = unname(files[idx[hit]])
  manifest
}

.feature_key = function(prefix, value) {
  value = tolower(.bundle_squish(value))
  value = gsub("[^a-z0-9]+", "_", value)
  value = gsub("^_+|_+$", "", value)
  value[value == ""] = "unknown"
  paste(prefix, value, sep = "__")
}

.feature_species_id = function(species) {
  out = tolower(.bundle_squish(species))
  out = gsub("[^a-z0-9]+", "_", out)
  out = gsub("^_+|_+$", "", out)
  paste0("species__", out)
}

.feature_species_universe = function(membership, species_universe = NULL) {
  supplied = if (is.data.frame(species_universe)) {
    names(species_universe) = .plant_normalize_column_names(
      names(species_universe)
    )
    if ("species" %in% names(species_universe)) {
      species_universe$species
    } else if (ncol(species_universe) > 0) {
      species_universe[[1]]
    } else {
      character()
    }
  } else {
    species_universe
  }
  supplied = unique(.uaf_non_empty(.bundle_squish(supplied)))
  observed = if (is.data.frame(membership) &&
                 "species" %in% names(membership)) {
    sort(unique(.uaf_non_empty(.bundle_squish(membership$species))))
  } else {
    character()
  }
  if (length(supplied) < 1) return(observed)
  c(supplied, setdiff(observed, supplied))
}

.feature_species_metadata = function(membership, species_universe = NULL,
                                      plant_metadata = NULL) {
  cols = c("species_id", "species", "accepted_species_name", "genus",
           "family", "taxonomy_family_status", "occurrence_count",
           "compound_count", "comparable_compound_count",
           "source_database_count", "direct_species_database_count",
           "direct_species_literature_count", "fallback_count",
           "candidate_only_count", "unresolved_or_review_required_count",
           "chemistry_record_status")
  species_universe = .feature_species_universe(membership, species_universe)
  if (length(species_universe) < 1) {
    return(.bundle_empty(cols))
  }
  if (!is.data.frame(membership)) membership = data.frame()
  if (!"species" %in% names(membership)) membership$species = character()
  grades = plantOccurrenceEvidenceGrade(membership)
  rows = lapply(species_universe, function(sp) {
    x = membership[membership$species == sp, , drop = FALSE]
    gx = grades[grades$species == sp, , drop = FALSE]
    mx = if (is.data.frame(plant_metadata) &&
             "species" %in% names(plant_metadata)) {
      plant_metadata[plant_metadata$species == sp, , drop = FALSE]
    } else {
      data.frame()
    }
    data.frame(
      species_id = .feature_species_id(sp),
      species = sp,
      accepted_species_name =
        .bundle_first_value(.bundle_first_non_empty(
          .bundle_col_or(x, "accepted_species_name", NA_character_),
          .bundle_col_or(mx, "accepted_species_name", NA_character_)
        )),
      genus = .bundle_first_value(.bundle_first_non_empty(
        .bundle_col_or(x, "genus", NA_character_),
        .bundle_col_or(mx, "genus", NA_character_), .bundle_genus(sp)
      )),
      family = .bundle_first_value(.bundle_first_non_empty(
        .bundle_col_or(x, "family", NA_character_),
        .bundle_col_or(mx, "family", NA_character_)
      )),
      taxonomy_family_status =
        .bundle_first_value(.bundle_first_non_empty(
          .bundle_col_or(x, "taxonomy_family_status", NA_character_),
          .bundle_col_or(mx, "family_status", NA_character_)
        )),
      occurrence_count = nrow(x),
      compound_count = length(unique(.uaf_non_empty(x$compound_id))),
      comparable_compound_count = length(unique(.uaf_non_empty(
        x$compound_id[.bundle_truthy(.bundle_col_or(x, "comparable_for_matrix",
                                                    NA_character_))]
      ))),
      source_database_count = length(unique(.uaf_non_empty(
        .bundle_col_or(x, "source_database", NA_character_)
      ))),
      direct_species_database_count = sum(
        gx$evidence_grade == "direct_species_database_record", na.rm = TRUE
      ),
      direct_species_literature_count = sum(
        gx$evidence_grade == "direct_species_literature_supported_record",
        na.rm = TRUE
      ),
      fallback_count = sum(
        gx$evidence_grade == "source_backed_genus_family_fallback",
        na.rm = TRUE
      ),
      candidate_only_count = sum(
        gx$evidence_grade == "pubtator_pubmed_candidate_only", na.rm = TRUE
      ),
      unresolved_or_review_required_count = sum(
        gx$evidence_grade == "unresolved_or_review_required", na.rm = TRUE
      ),
      chemistry_record_status = if (nrow(x) > 0) {
        "records_present"
      } else {
        "no_records_in_membership"
      },
      stringsAsFactors = FALSE
    )
  })
  out = .bundle_bind(rows)
  out[, cols, drop = FALSE]
}

.project_read_result = function(plant_result) {
  if (is.null(plant_result)) return(NULL)
  if (is.character(plant_result) && length(plant_result) == 1 &&
      file.exists(plant_result)) {
    return(readRDS(plant_result))
  }
  plant_result
}

.project_membership_from_inputs = function(result = NULL,
                                           plant_compounds = NULL) {
  if (!is.null(plant_compounds)) {
    x = if (is.character(plant_compounds) && length(plant_compounds) == 1 &&
        file.exists(plant_compounds)) {
      utils::read.csv(plant_compounds, stringsAsFactors = FALSE,
                      check.names = FALSE)
    } else {
      plant_compounds
    }
    if (!is.data.frame(x)) {
      stop("`plant_compounds` must be a data frame or CSV path.",
           call. = FALSE)
    }
    if (all(c("species", "compound_name") %in%
            .plant_normalize_column_names(names(x)))) {
      raw = as.data.frame(x, stringsAsFactors = FALSE)
      names(raw) = .plant_normalize_column_names(names(raw))
      out = standardizePlantCompoundIntake(raw)
      if (nrow(out) == nrow(raw)) {
        for (col in setdiff(names(raw), names(out))) out[[col]] = raw[[col]]
      }
      return(out)
    }
    return(as.data.frame(x, stringsAsFactors = FALSE))
  }
  if (is.list(result) && is.data.frame(result$PlantCompoundOccurrences)) {
    return(as.data.frame(result$PlantCompoundOccurrences,
                         stringsAsFactors = FALSE))
  }
  if (is.list(result) && is.data.frame(result$PlantCompoundMembership)) {
    return(as.data.frame(result$PlantCompoundMembership,
                         stringsAsFactors = FALSE))
  }
  NULL
}

.project_write_minimal_bundle = function(bundle_dir, membership, plants,
                                         metadata = NULL, project_id = NULL) {
  if (dir.exists(bundle_dir) || file.exists(bundle_dir)) {
    unlink(bundle_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(bundle_dir, recursive = TRUE, showWarnings = FALSE)
  .categorate_write_csv_file(
    data.frame(
      Table = "PlantCompoundMembership",
      SheetName = "PlantCompoundMembership",
      FileName = "03_PlantCompoundMembership.csv",
      RowCount = nrow(membership),
      ColumnCount = ncol(membership),
      Format = "csv",
      OutputPath = ".",
      PathType = "bundle_relative",
      CreatedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      ProjectID = .uaf_first_non_empty_text(project_id, NA_character_),
      stringsAsFactors = FALSE
    ),
    file.path(bundle_dir, "01_ExportManifest.csv")
  )
  .categorate_write_csv_file(membership,
                             file.path(bundle_dir,
                                       "03_PlantCompoundMembership.csv"))
  if (is.data.frame(plants) && nrow(plants) > 0) {
    .categorate_write_csv_file(plants,
                               file.path(bundle_dir,
                                         "07_ProjectPlantList.csv"))
  }
  if (is.data.frame(metadata) && nrow(metadata) > 0) {
    .categorate_write_csv_file(metadata,
                               file.path(bundle_dir,
                                         "07b_ProjectPlantMetadata.csv"))
  }
  .bundle_refresh_manifest(bundle_dir, project_id = project_id)
  invisible(bundle_dir)
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
    out$identity_issue_type = .bundle_first_non_empty(
      .bundle_col_or(out, "identity_issue_type", NA_character_),
      .bundle_index(resolved, idx, "identity_issue_type")
    )
    out$identity_review_reason = .bundle_first_non_empty(
      .bundle_col_or(out, "identity_review_reason", NA_character_),
      .bundle_index(resolved, idx, "review_reason")
    )
    out$identity_recommended_action = .bundle_first_non_empty(
      .bundle_col_or(out, "identity_recommended_action", NA_character_),
      .bundle_index(resolved, idx, "recommended_action")
    )
    out$identity_ambiguity_flag = .bundle_first_non_empty(
      .bundle_col_or(out, "identity_ambiguity_flag", NA_character_),
      .bundle_index(resolved, idx, "ambiguity_flag")
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

  plant_part_known = .bundle_known(
    .bundle_col_or(out, "plant_part_group", NA_character_)
  )
  tissue_known = .bundle_known(
    .bundle_col_or(out, "tissue_group", NA_character_)
  )
  analytical_method_known = .bundle_analytical_method_known(
    .bundle_col_or(out, "method_group", NA_character_)
  )
  source_provenance_record = .bundle_source_provenance_record(
    .bundle_col_or(out, "method_group", NA_character_)
  )
  biological_context_known = plant_part_known | tissue_known
  out$context_known_record = .bundle_yes_no(biological_context_known)
  out$biological_context_known = .bundle_yes_no(biological_context_known)
  out$plant_part_known = .bundle_yes_no(plant_part_known)
  out$tissue_known = .bundle_yes_no(tissue_known)
  out$method_known = .bundle_yes_no(analytical_method_known)
  out$analytical_method_known = .bundle_yes_no(analytical_method_known)
  out$source_provenance_record = .bundle_yes_no(source_provenance_record)

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
    "identity_issue_type", "identity_review_reason",
    "identity_recommended_action", "identity_ambiguity_flag",
    "context_known_record", "biological_context_known",
    "plant_part_known", "tissue_known", "method_known",
    "analytical_method_known", "source_provenance_record",
    "has_fingerprint", "metabolism_domain",
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
    if (!col %in% names(x)) x[[col]] = rep(NA, nrow(x))
  }
  if (nrow(x) < 1) {
    for (col in c("unordered_pair_key", "support_tier",
                  "low_support_caution", "support_note")) {
      if (!col %in% names(x)) x[[col]] = character()
    }
    return(x)
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
           "tissue_known_fraction", "analytical_method_known_fraction",
           "biological_context_known_fraction",
           "context_known_record_fraction",
           "source_provenance_record_fraction",
           "identity_review_compound_fraction",
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
      tissue_known_fraction = .bundle_fraction(
        .bundle_col_or(x, "tissue_known", NA_character_)
      ),
      analytical_method_known_fraction = .bundle_fraction(
        .bundle_col_or(x, "analytical_method_known", NA_character_)
      ),
      biological_context_known_fraction = .bundle_fraction(
        .bundle_col_or(x, "biological_context_known", NA_character_)
      ),
      context_known_record_fraction = .bundle_fraction(
        .bundle_col_or(x, "context_known_record", NA_character_)
      ),
      source_provenance_record_fraction = .bundle_fraction(
        .bundle_col_or(x, "source_provenance_record", NA_character_)
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
  schema = uafRSchemaMetadata(
    workflow_name = "finalizePlantChemistryAnalysisBundle",
    workflow_parameters = list(project_id = project_id)
  )
  rows = lapply(csv_files, function(file_name) {
    table = .bundle_table_from_file(file_name, existing)
    file = file.path(path, file_name)
    dat = tryCatch(utils::read.csv(file,
                                   stringsAsFactors = FALSE,
                                   check.names = FALSE),
                   error = function(e) data.frame())
    artifact = .bundle_artifact_metadata(file, file_name)
    data.frame(
      Table = table,
      SheetName = substr(gsub("[^A-Za-z0-9_]+", "_", table), 1, 31),
      FileName = file_name,
      RowCount = nrow(dat),
      ColumnCount = ncol(dat),
      Format = "csv",
      OutputPath = ".",
      PathType = "bundle_relative",
      ArtifactPath = file_name,
      FileSizeBytes = artifact$size,
      ChecksumAlgorithm = artifact$algorithm,
      ArtifactChecksum = artifact$checksum,
      ChecksumStatus = artifact$status,
      CreatedAt = created_at,
      ProjectID = .uaf_first_non_empty_text(project_id, NA_character_),
      uafR_schema_version = schema$uafR_schema_version[[1]],
      uafR_package_version = schema$uafR_package_version[[1]],
      created_at = schema$created_at[[1]],
      workflow_name = schema$workflow_name[[1]],
      workflow_parameters = schema$workflow_parameters[[1]],
      stringsAsFactors = FALSE
    )
  })
  manifest = .bundle_bind(rows)
  markdown = c("README.md", "METHODS_TEXT.md")
  md_rows = lapply(markdown[file.exists(file.path(path, markdown))],
                   function(file_name) {
    artifact = .bundle_artifact_metadata(file.path(path, file_name),
                                         file_name)
    data.frame(
      Table = sub("[.].*$", "", file_name),
      SheetName = sub("[.].*$", "", file_name),
      FileName = file_name,
      RowCount = NA_integer_,
      ColumnCount = NA_integer_,
      Format = "markdown",
      OutputPath = ".",
      PathType = "bundle_relative",
      ArtifactPath = file_name,
      FileSizeBytes = artifact$size,
      ChecksumAlgorithm = artifact$algorithm,
      ArtifactChecksum = artifact$checksum,
      ChecksumStatus = artifact$status,
      CreatedAt = created_at,
      ProjectID = .uaf_first_non_empty_text(project_id, NA_character_),
      uafR_schema_version = schema$uafR_schema_version[[1]],
      uafR_package_version = schema$uafR_package_version[[1]],
      created_at = schema$created_at[[1]],
      workflow_name = schema$workflow_name[[1]],
      workflow_parameters = schema$workflow_parameters[[1]],
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

.bundle_artifact_metadata = function(file, file_name = basename(file)) {
  is_manifest = identical(file_name, "01_ExportManifest.csv")
  exists = file.exists(file)
  size = if (exists && !isTRUE(file.info(file)$isdir)) {
    as.numeric(file.info(file)$size)
  } else {
    NA_real_
  }
  if (is_manifest) {
    return(list(size = NA_real_, algorithm = NA_character_,
                checksum = NA_character_, status = "not_self_hashed"))
  }
  if (!exists) {
    return(list(size = NA_real_, algorithm = "MD5",
                checksum = NA_character_, status = "missing"))
  }
  checksum = tryCatch(unname(tools::md5sum(file)[[1]]),
                      error = function(error) NA_character_)
  list(
    size = size,
    algorithm = "MD5",
    checksum = checksum,
    status = if (.bundle_known(checksum)) "computed" else "error"
  )
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
    "15_PlantChemistryMissingSpecies.csv" = "PlantChemistryMissingSpecies",
    "16_EvidenceGradeSummary.csv" = "EvidenceGradeSummary",
    "17_ReviewRequiredOccurrences.csv" = "ReviewRequiredOccurrences",
    "18_FeatureComparisonGroupCountMatrix.csv" =
      "FeatureComparisonGroupCountMatrix",
    "19_FeatureComparisonScopeCountMatrix.csv" =
      "FeatureComparisonScopeCountMatrix",
    "20_FeatureSourceCoverageMatrix.csv" = "FeatureSourceCoverageMatrix",
    "21_FeatureSpeciesMetadata.csv" = "FeatureSpeciesMetadata",
    "22_ComparableScopeTanimotoSummary.csv" =
      "ComparableScopeTanimotoSummary",
    "23_ComparableGroupTanimotoSummary.csv" =
      "ComparableGroupTanimotoSummary",
    "24_FeatureEvidenceGradeCountMatrix.csv" =
      "FeatureEvidenceGradeCountMatrix",
    "25_FeaturePlantPartCountMatrix.csv" =
      "FeaturePlantPartCountMatrix",
    "26_FeatureTissueCountMatrix.csv" =
      "FeatureTissueCountMatrix",
    "27_FeatureMethodCountMatrix.csv" =
      "FeatureMethodCountMatrix",
    "28_FeatureMatrixManifest.csv" = "FeatureMatrixManifest",
    "07_ProjectPlantList.csv" = "ProjectPlantList",
    "07b_ProjectPlantMetadata.csv" = "ProjectPlantMetadata"
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
    "- `01_ExportManifest.csv`: bundle-relative file names, row/column counts, sizes, MD5 checksums, schema/package versions, and creation metadata. Its own checksum is intentionally omitted to avoid a circular hash.",
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
    ExportManifest = c("Table", "FileName", "RowCount", "ColumnCount",
                       "OutputPath", "PathType", "ArtifactPath",
                       "FileSizeBytes", "ChecksumAlgorithm",
                       "ArtifactChecksum", "ChecksumStatus",
                       "uafR_schema_version", "uafR_package_version",
                       "workflow_name"),
    PlantCompoundMembership = c("species", "compound_id", "compound_name",
                                "source_database", "evidence_tier"),
    PlantCompoundMembershipEnriched = c(
      "species", "compound_id", "compound_name", "comparison_scope",
      "comparison_group", "has_fingerprint", "biological_context_known",
      "analytical_method_known", "source_provenance_record"
    ),
    PlantPairTanimotoSummary = c("species_a", "species_b",
                                 "compound_pair_count", "mean_tanimoto",
                                 "unordered_pair_key", "support_tier"),
    PlantChemistrySummary = c("species", "reported_compound_count",
                              "fingerprinted_compound_count",
                              "biological_context_known_fraction",
                              "chemistry_data_quality_tier"),
    EvidenceGradeSummary = c("evidence_grade", "occurrence_count",
                             "review_required_count"),
    ReviewRequiredOccurrences = c("species", "compound_id",
                                  "evidence_grade", "review_category",
                                  "review_reason", "recommended_action"),
    FeatureComparisonGroupCountMatrix = c("species_id", "species"),
    FeatureComparisonScopeCountMatrix = c("species_id", "species"),
    FeatureSourceCoverageMatrix = c("species_id", "species"),
    FeatureSpeciesMetadata = c("species_id", "species", "compound_count",
                               "chemistry_record_status"),
    FeatureEvidenceGradeCountMatrix = c("species_id", "species"),
    FeaturePlantPartCountMatrix = c("species_id", "species"),
    FeatureTissueCountMatrix = c("species_id", "species"),
    FeatureMethodCountMatrix = c("species_id", "species"),
    FeatureMatrixManifest = c("Table", "FileName", "RowCount",
                              "ColumnCount", "FeatureField", "Mode",
                              "SpeciesUniverseCount"),
    ComparableScopeTanimotoSummary = c("species_a", "species_b",
                                       "comparison_scope",
                                       "compound_pair_count",
                                       "mean_tanimoto", "support_tier"),
    ComparableGroupTanimotoSummary = c("species_a", "species_b",
                                       "comparison_group",
                                       "compound_pair_count",
                                       "mean_tanimoto", "support_tier")
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

.bundle_validate_export_artifacts = function(path, manifest) {
  cols = c("Table", "FileName", "ArtifactPath", "PathType", "FileExists",
           "ExpectedSizeBytes", "ActualSizeBytes", "SizeMatches",
           "ChecksumAlgorithm", "ExpectedChecksum", "ActualChecksum",
           "ChecksumMatches", "PortablePath", "Status", "Message")
  if (!is.data.frame(manifest) || nrow(manifest) < 1) {
    return(data.frame(
      Table = NA_character_, FileName = NA_character_,
      ArtifactPath = NA_character_, PathType = NA_character_,
      FileExists = "No", ExpectedSizeBytes = NA_real_,
      ActualSizeBytes = NA_real_, SizeMatches = "No",
      ChecksumAlgorithm = NA_character_, ExpectedChecksum = NA_character_,
      ActualChecksum = NA_character_, ChecksumMatches = "No",
      PortablePath = "No", Status = "fail",
      Message = "Bundle export manifest is missing or empty.",
      stringsAsFactors = FALSE
    )[, cols, drop = FALSE])
  }
  required = c("Table", "FileName", "OutputPath", "PathType",
               "ArtifactPath", "FileSizeBytes", "ChecksumAlgorithm",
               "ArtifactChecksum", "ChecksumStatus")
  missing = setdiff(required, names(manifest))
  if (length(missing) > 0) {
    return(data.frame(
      Table = NA_character_, FileName = NA_character_,
      ArtifactPath = NA_character_, PathType = NA_character_,
      FileExists = "No", ExpectedSizeBytes = NA_real_,
      ActualSizeBytes = NA_real_, SizeMatches = "No",
      ChecksumAlgorithm = NA_character_, ExpectedChecksum = NA_character_,
      ActualChecksum = NA_character_, ChecksumMatches = "No",
      PortablePath = "No", Status = "warn",
      Message = paste("Legacy manifest lacks portable artifact fields:",
                      paste(missing, collapse = ", ")),
      stringsAsFactors = FALSE
    )[, cols, drop = FALSE])
  }
  rows = lapply(seq_len(nrow(manifest)), function(i) {
    row = manifest[i, , drop = FALSE]
    file_name = .uaf_first_non_empty_text(row$FileName, NA_character_)
    artifact_path = .uaf_first_non_empty_text(row$ArtifactPath, file_name)
    output_path = .uaf_first_non_empty_text(row$OutputPath, NA_character_)
    path_type = .uaf_first_non_empty_text(row$PathType, NA_character_)
    portable = identical(output_path, ".") &&
      identical(path_type, "bundle_relative") &&
      .bundle_safe_relative_artifact_path(artifact_path) &&
      identical(artifact_path, file_name)
    file = if (.bundle_safe_relative_artifact_path(artifact_path)) {
      file.path(path, artifact_path)
    } else {
      NA_character_
    }
    exists = length(file) == 1 && !is.na(file) && file.exists(file)
    actual_size = if (exists && !isTRUE(file.info(file)$isdir)) {
      as.numeric(file.info(file)$size)
    } else {
      NA_real_
    }
    expected_size = suppressWarnings(as.numeric(row$FileSizeBytes))
    self_manifest = identical(file_name, "01_ExportManifest.csv")
    size_matches = if (self_manifest) {
      TRUE
    } else {
      exists && !is.na(expected_size) && identical(actual_size, expected_size)
    }
    expected_checksum = .uaf_first_non_empty_text(row$ArtifactChecksum,
                                                   NA_character_)
    checksum_algorithm = .uaf_first_non_empty_text(row$ChecksumAlgorithm,
                                                   NA_character_)
    checksum_status = .uaf_first_non_empty_text(row$ChecksumStatus,
                                                NA_character_)
    actual_checksum = if (exists && !self_manifest &&
                           identical(checksum_status, "computed")) {
      tryCatch(unname(tools::md5sum(file)[[1]]),
               error = function(error) NA_character_)
    } else {
      NA_character_
    }
    algorithm_valid = if (self_manifest) {
      is.na(checksum_algorithm) || !nzchar(checksum_algorithm)
    } else {
      identical(toupper(checksum_algorithm), "MD5")
    }
    checksum_matches = if (self_manifest && algorithm_valid &&
                            identical(checksum_status, "not_self_hashed")) {
      TRUE
    } else {
      exists && algorithm_valid && identical(checksum_status, "computed") &&
        .bundle_known(expected_checksum) &&
        identical(actual_checksum, expected_checksum)
    }
    messages = character()
    if (!portable) messages = c(messages, "artifact path is not portable")
    if (!exists) messages = c(messages, "artifact file is missing")
    if (!size_matches) messages = c(messages, "artifact size differs")
    if (!algorithm_valid) messages = c(messages, "checksum algorithm is invalid")
    if (!checksum_matches) messages = c(messages, "artifact checksum differs")
    ok = portable && exists && size_matches && checksum_matches
    data.frame(
      Table = as.character(row$Table), FileName = file_name,
      ArtifactPath = artifact_path, PathType = path_type,
      FileExists = .bundle_yes_no(exists),
      ExpectedSizeBytes = expected_size, ActualSizeBytes = actual_size,
      SizeMatches = .bundle_yes_no(size_matches),
      ChecksumAlgorithm = checksum_algorithm,
      ExpectedChecksum = expected_checksum,
      ActualChecksum = actual_checksum,
      ChecksumMatches = .bundle_yes_no(checksum_matches),
      PortablePath = .bundle_yes_no(portable),
      Status = if (ok) "pass" else "fail",
      Message = if (ok) {
        if (self_manifest) {
          "Manifest path is portable; self-checksum is intentionally omitted."
        } else {
          "Artifact path, size, and checksum are valid."
        }
      } else {
        paste(unique(messages), collapse = "; ")
      },
      stringsAsFactors = FALSE
    )
  })
  .bundle_bind(rows)[, cols, drop = FALSE]
}

.bundle_safe_relative_artifact_path = function(x) {
  x = as.character(x)
  length(x) == 1 && !is.na(x) && nzchar(x) &&
    !grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\|~[/\\\\])", x) &&
    !grepl("(^|[/\\\\])[.][.]($|[/\\\\])", x) &&
    identical(basename(x), x)
}

.bundle_validate_manifest_references = function(path, manifest) {
  cols = c("ManifestTable", "ReferencedTable", "FileName", "FileExists",
           "ListedInExportManifest", "ExpectedRows", "ActualRows",
           "ExpectedColumns", "ActualColumns", "SpeciesUniverseMatches",
           "Status", "Message")
  if (!is.data.frame(manifest) || nrow(manifest) < 1 ||
      !"FeatureMatrixManifest" %in% manifest$Table) {
    return(.bundle_empty(cols))
  }
  nested = .bundle_read_table(path, manifest, "FeatureMatrixManifest")
  required = c("Table", "FileName", "RowCount", "ColumnCount")
  missing = setdiff(required, names(nested))
  if (length(missing) > 0) {
    return(data.frame(
      ManifestTable = "FeatureMatrixManifest",
      ReferencedTable = NA_character_, FileName = NA_character_,
      FileExists = "No", ListedInExportManifest = "No",
      ExpectedRows = NA_integer_, ActualRows = NA_integer_,
      ExpectedColumns = NA_integer_, ActualColumns = NA_integer_,
      SpeciesUniverseMatches = "No", Status = "fail",
      Message = paste("Missing nested manifest columns:",
                      paste(missing, collapse = ", ")),
      stringsAsFactors = FALSE
    ))
  }
  metadata_row = which(nested$Table == "SpeciesMetadata")
  expected_species = character()
  if (length(metadata_row) > 0) {
    metadata_file = file.path(path, nested$FileName[[metadata_row[[1]]]])
    metadata = tryCatch(
      utils::read.csv(metadata_file, stringsAsFactors = FALSE,
                      check.names = FALSE),
      error = function(error) data.frame()
    )
    if ("species_id" %in% names(metadata)) {
      expected_species = as.character(metadata$species_id)
    }
  }
  rows = lapply(seq_len(nrow(nested)), function(i) {
    file_name = as.character(nested$FileName[[i]])
    file = file.path(path, file_name)
    exists = file.exists(file)
    dat = if (exists) {
      tryCatch(utils::read.csv(file, stringsAsFactors = FALSE,
                               check.names = FALSE),
               error = function(error) data.frame())
    } else {
      data.frame()
    }
    listed = "FileName" %in% names(manifest) && file_name %in%
      manifest$FileName
    expected_rows = suppressWarnings(as.integer(nested$RowCount[[i]]))
    expected_cols = suppressWarnings(as.integer(nested$ColumnCount[[i]]))
    actual_rows = if (exists) nrow(dat) else NA_integer_
    actual_cols = if (exists) ncol(dat) else NA_integer_
    species_match = exists && "species_id" %in% names(dat) &&
      length(expected_species) > 0 &&
      identical(as.character(dat$species_id), expected_species)
    ok = exists && listed && !is.na(expected_rows) &&
      identical(actual_rows, expected_rows) && !is.na(expected_cols) &&
      identical(actual_cols, expected_cols) && species_match
    messages = character()
    if (!exists) messages = c(messages, "referenced file is missing")
    if (!listed) messages = c(messages, "file is absent from export manifest")
    if (exists && !identical(actual_rows, expected_rows)) {
      messages = c(messages, "row count differs from nested manifest")
    }
    if (exists && !identical(actual_cols, expected_cols)) {
      messages = c(messages, "column count differs from nested manifest")
    }
    if (exists && !species_match) {
      messages = c(messages, "species universe/order differs from SpeciesMetadata")
    }
    data.frame(
      ManifestTable = "FeatureMatrixManifest",
      ReferencedTable = as.character(nested$Table[[i]]),
      FileName = file_name,
      FileExists = .bundle_yes_no(exists),
      ListedInExportManifest = .bundle_yes_no(listed),
      ExpectedRows = expected_rows, ActualRows = actual_rows,
      ExpectedColumns = expected_cols, ActualColumns = actual_cols,
      SpeciesUniverseMatches = .bundle_yes_no(species_match),
      Status = if (ok) "pass" else "fail",
      Message = if (length(messages) > 0) {
        paste(messages, collapse = "; ")
      } else {
        "Referenced feature artifact is present and internally consistent."
      },
      stringsAsFactors = FALSE
    )
  })
  .bundle_bind(rows)[, cols, drop = FALSE]
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

.bundle_analytical_method_known = function(x) {
  x = tolower(.bundle_squish(as.character(x)))
  .bundle_known(x) & !x %in% c("database_record", "literature_curation")
}

.bundle_source_provenance_record = function(x) {
  x = tolower(.bundle_squish(as.character(x)))
  !is.na(x) & x %in% c("database_record", "literature_curation")
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
         EvidenceGradeSummary = "Counts plant-compound occurrence rows by conservative evidence grade.",
         ReviewRequiredOccurrences = "Occurrence rows requiring curation before strong downstream use.",
         FeatureComparisonGroupCountMatrix = "Species by comparable chemistry group count matrix.",
         FeatureComparisonScopeCountMatrix = "Species by comparable chemistry scope count matrix.",
         FeatureSourceCoverageMatrix = "Species by source database count matrix.",
         FeatureSpeciesMetadata = "Species-level quality, coverage, taxonomy, and evidence metadata.",
         FeatureEvidenceGradeCountMatrix = "Species by conservative evidence-grade count matrix.",
         FeaturePlantPartCountMatrix = "Species by source-backed plant-part context count matrix.",
         FeatureTissueCountMatrix = "Species by source-backed tissue context count matrix.",
         FeatureMethodCountMatrix = "Species by source-backed method context count matrix.",
         FeatureMatrixManifest = "Lists model-ready feature matrices, modes, filters, and caveats.",
         ComparableScopeTanimotoSummary = "Plant-pair Tanimoto summary filtered to matching comparable chemistry scopes.",
         ComparableGroupTanimotoSummary = "Plant-pair Tanimoto summary filtered to matching comparable chemistry groups.",
         ProjectPlantList = "Plant list supplied to a project runner.",
         ProjectPlantMetadata = "Standardized project metadata supplied to a project runner.",
         SourceCoverageSummary = "Source-level coverage summary from categorate batches.",
         ValidationOverview = "Bundle-level validation and export readiness summary.",
         "uafR plant chemistry bundle table.")
}

.bundle_column_required = function(table, col) {
  req = list(
    ExportManifest = c("Table", "FileName", "RowCount", "ColumnCount",
                       "OutputPath", "PathType", "ArtifactPath",
                       "FileSizeBytes", "ChecksumAlgorithm",
                       "ArtifactChecksum", "ChecksumStatus",
                       "uafR_schema_version", "uafR_package_version",
                       "workflow_name"),
    PlantCompoundMembership = c("species", "compound_id", "compound_name"),
    PlantCompoundMembershipEnriched = c("species", "compound_id",
                                        "compound_name", "comparison_scope",
                                        "biological_context_known"),
    PlantPairTanimotoSummary = c("species_a", "species_b",
                                 "unordered_pair_key"),
    PlantChemistrySummary = c("species", "reported_compound_count",
                              "biological_context_known_fraction"),
    EvidenceGradeSummary = c("evidence_grade", "occurrence_count"),
    ReviewRequiredOccurrences = c("species", "compound_id",
                                  "evidence_grade", "review_category",
                                  "recommended_action"),
    FeatureComparisonGroupCountMatrix = c("species_id", "species"),
    FeatureComparisonScopeCountMatrix = c("species_id", "species"),
    FeatureSourceCoverageMatrix = c("species_id", "species"),
    FeatureSpeciesMetadata = c("species_id", "species", "compound_count",
                               "chemistry_record_status"),
    FeatureEvidenceGradeCountMatrix = c("species_id", "species"),
    FeaturePlantPartCountMatrix = c("species_id", "species"),
    FeatureTissueCountMatrix = c("species_id", "species"),
    FeatureMethodCountMatrix = c("species_id", "species"),
    FeatureMatrixManifest = c("Table", "FileName", "RowCount",
                              "ColumnCount", "FeatureField", "Mode",
                              "SpeciesUniverseCount"),
    ComparableScopeTanimotoSummary = c("species_a", "species_b",
                                       "comparison_scope",
                                       "compound_pair_count"),
    ComparableGroupTanimotoSummary = c("species_a", "species_b",
                                       "comparison_group",
                                       "compound_pair_count")
  )
  table %in% names(req) && col %in% req[[table]]
}

.bundle_allowed_values = function(col) {
  switch(col,
         comparable_for_matrix = "Yes; No",
         has_fingerprint = "Yes; No",
         context_known_record = "Yes; No",
         biological_context_known = "Yes; No",
         plant_part_known = "Yes; No",
         tissue_known = "Yes; No",
         method_known = "Yes; No",
         analytical_method_known = "Yes; No",
         source_provenance_record = "Yes; No",
         matched_in_lotus = "Yes; No",
         source_backed = "Yes; No",
         structure_resolved = "Yes; No",
         comparable_for_analysis = "Yes; No",
         review_required = "Yes; No",
         evidence_review_required = "Yes; No",
         identity_review_required = "Yes; No",
         structure_review_required = "Yes; No",
         context_review_required = "Yes; No",
         comparability_review_required = "Yes; No",
         citation_review_required = "Yes; No",
         low_support_caution = "Yes; No",
         accepted_name_status = "accepted_name_supplied; accepted_name_not_supplied",
         family_status = "family_supplied; family_not_supplied",
         evidence_grade = "direct_species_database_record; direct_species_literature_supported_record; source_backed_genus_family_fallback; pubtator_pubmed_candidate_only; unresolved_or_review_required; excluded_by_review",
         support_tier = "very_low; low; moderate; high; unknown",
         ExportReadyStatus = "pass; warn; fail",
         chemistry_record_status = "records_present; no_records_in_membership",
         "")
}

.bundle_column_note = function(col) {
  switch(col,
         ArtifactPath = "Bundle-relative artifact name; resolve against the bundle directory.",
         FileSizeBytes = "Artifact size recorded when the manifest was refreshed.",
         ArtifactChecksum = "MD5 integrity value for the artifact; the export manifest omits its own checksum to avoid a circular hash.",
         ChecksumStatus = "Reports whether a checksum was computed or intentionally omitted for the manifest itself.",
         PathType = "Distinguishes portable bundle-relative paths from external reference basenames.",
         comparison_scope = "Use to keep downstream comparisons biologically comparable.",
         comparison_group = "Higher-level chemistry grouping within a comparison scope.",
         evidence_tier = "Source evidence level; literature co-mentions should not be treated as confirmed occurrence unless curated.",
         evidence_grade = "Conservative uafR analysis tier derived from source, rank, evidence tier, and optional review decisions.",
         recommended_use = "Suggested downstream use based on evidence grade, structure resolution, and comparability.",
         review_category = "Semicolon-delimited evidence, identity, structure, biological-context, comparability, or citation review categories.",
         recommended_action = "Specific curation or exclusion action required before the affected downstream analysis.",
         analysis_filter_key = "Convenience key for filtering source-backed, structure-resolved, comparable rows.",
         species_id = "Stable sanitized species identifier for joins and model-ready matrices.",
         biological_context_known = "Yes only when a source reports plant-part or tissue context; database provenance alone does not qualify.",
         analytical_method_known = "Yes only for a reported analytical method; generic database or literature provenance does not qualify.",
         source_provenance_record = "Identifies generic database or literature provenance separately from biological and analytical context.",
         SpeciesUniverseCount = "Complete species count that every feature matrix is required to retain in the same order.",
         comparison_value = "Scope or group value used for a filtered comparable chemistry summary.",
         comparison_filter = "Human-readable filter used to create the filtered comparison table.",
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
         evidence_grade = "Evidence grades are conservative analysis tiers, not experimental confirmation.",
         source_backed = "Source-backed means a public/source record exists; it does not prove occurrence in project samples.",
         comparable_for_analysis = "Comparable flags are derived from available classification evidence and should be reviewed for the project question.",
         comparison_value = "Filtered Tanimoto summaries are only as strong as the compound classifications and support counts behind them.",
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
