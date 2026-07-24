#' Pairwise PubChem fingerprint Tanimoto similarity for chemicals
#'
#' @description
#' `chemicalTanimotoSimilarity()` resolves chemicals to PubChem CIDs, fetches
#' PubChem Fingerprint2D bit strings, decodes them with `ChemmineR`, and
#' computes pairwise Tanimoto similarity. It accepts a plain character vector or
#' a data frame with compound identifiers. When a grouping column is supplied,
#' it also returns group-level summaries that are suitable for joining to other
#' sample, species, treatment, or site similarity tables.
#'
#' @param compounds Character vector or data frame. Character vectors are
#' treated as compound names. Data frames can contain compound names, PubChem
#' CIDs, InChIKeys, SMILES, and optional grouping columns.
#' @param compound_col Column containing compound names when `compounds` is a
#' data frame. If `NULL`, common names such as `compound_name`, `Query`, or
#' `Chemical` are detected.
#' @param cid_col Column containing PubChem CIDs. If `NULL`, common names such
#' as `CID`, `pubchem_cid`, or `cid` are detected.
#' @param inchikey_col Column containing InChIKeys. If `NULL`, common names such
#' as `InChIKey` or `inchikey` are detected.
#' @param smiles_col Optional SMILES column retained as metadata. PubChem
#' Fingerprint2D values are still fetched from PubChem so the Tanimoto
#' definition is consistent across workflows.
#' @param compound_id_col Optional stable compound identifier column. If absent,
#' identifiers are derived from InChIKey, PubChem CID, or compound name.
#' @param group_cols Optional grouping columns. For plants, use `"species"` or
#' call `plantChemicalTanimotoSimilarity()`.
#' @param thresholds Numeric Tanimoto thresholds to count in group summaries.
#' @param top_n_pairs Number of top compound-pair explanations to include in
#' group summaries.
#' @param return_compound_pairs Logical. If `TRUE`, include or write the
#' compound-compound pair table.
#' @param return_group_compound_pairs Logical. If `TRUE` and `group_cols` are
#' supplied, include or write all cross-group compound-pair rows.
#' @param out_dir Optional output directory. When supplied, large pairwise
#' tables are written as compressed CSV files and the result records their
#' paths.
#' @param max_in_memory_pairs Maximum pair rows allowed in memory when
#' `out_dir` is `NULL`.
#' @param pair_block_size Number of focal rows per pairwise computation block.
#' @param pair_shard_rows Maximum rows per compressed pairwise CSV shard.
#' Use `Inf` (the default) to preserve one file per requested pair table.
#' @param cache Logical. If `TRUE`, PubChem JSON responses are cached.
#' @param cache_dir Cache directory. Defaults to the uafR user cache.
#' @param throttle Seconds to wait after uncached PubChem requests.
#' @param refresh Logical. If `TRUE`, ignore cached PubChem responses.
#' @param request_fun Optional request function for tests. It receives a URL and
#' returns parsed JSON or JSON text.
#' @param name_fallback Logical. If `TRUE`, unresolved compounds may be queried
#' by compound name after CID and InChIKey attempts fail.
#' @param progress_fun Optional callback receiving named progress-event lists.
#' It is intended for durable project runners and is ignored by default.
#' @param progress_every Emit identity and pair-block progress after this many
#' completed items or blocks.
#' @param service_busy_limit Number of consecutive PubChem HTTP 429/503
#' responses allowed before a `uaf_pubchem_service_busy` condition stops the
#' run. The default `Inf` preserves ordinary interactive behavior.
#'
#' @return A list with class `"uaf_tanimoto_similarity"` containing
#' `CompoundInput`, `CompoundResolution`, `PubChemFingerprints`,
#' `CompoundIdentityMap`, `ExcludedCompoundIdentities`, `CanonicalCompounds`,
#' `GroupCompoundMembershipEvidence`, `CompoundTanimoto`,
#' `GroupCompoundMembership`, `GroupPairTanimotoSummary`,
#' `GroupCompoundTanimoto`, `ComparableScopeTanimotoSummary`,
#' `ComparableGroupTanimotoSummary`, `ExportManifest`, and `Provenance`. The
#' identity tables retain every input mapping while analysis membership is
#' collapsed to verified canonical PubChem structures. Comparable summaries
#' require source-backed comparison labels and exclude unknown or explicitly
#' non-comparable membership by default.
#'
#' @examples
#' \dontrun{
#' compounds = data.frame(
#'   species = c("Plant A", "Plant A", "Plant B"),
#'   compound_name = c("caffeine", "theobromine", "aspirin")
#' )
#' sim = chemicalTanimotoSimilarity(compounds, group_cols = "species")
#' sim$GroupPairTanimotoSummary
#' }
#'
#' @export
chemicalTanimotoSimilarity = function(compounds,
                                      compound_col = NULL,
                                      cid_col = NULL,
                                      inchikey_col = NULL,
                                      smiles_col = NULL,
                                      compound_id_col = NULL,
                                      group_cols = NULL,
                                      thresholds = c(0.5, 0.7, 0.85, 0.95),
                                      top_n_pairs = 5,
                                      return_compound_pairs = TRUE,
                                      return_group_compound_pairs = FALSE,
                                      out_dir = NULL,
                                      max_in_memory_pairs = 1000000,
                                      pair_block_size = 250,
                                      pair_shard_rows = Inf,
                                      cache = TRUE,
                                      cache_dir = NULL,
                                      throttle = 0.2,
                                      refresh = FALSE,
                                      request_fun = NULL,
                                      name_fallback = TRUE,
                                      progress_fun = NULL,
                                      progress_every = 25,
                                      service_busy_limit = Inf) {
  thresholds = suppressWarnings(as.numeric(thresholds))
  thresholds = thresholds[!is.na(thresholds)]
  if (length(thresholds) < 1) thresholds = c(0.5, 0.7, 0.85, 0.95)
  pair_block_size = suppressWarnings(as.integer(pair_block_size))
  if (length(pair_block_size) != 1L || is.na(pair_block_size) ||
      pair_block_size < 1L) {
    stop("`pair_block_size` must be a positive integer.", call. = FALSE)
  }
  pair_shard_rows = suppressWarnings(as.numeric(pair_shard_rows))
  if (length(pair_shard_rows) != 1L || is.na(pair_shard_rows) ||
      pair_shard_rows < 1) {
    stop("`pair_shard_rows` must be a positive number or Inf.",
         call. = FALSE)
  }
  progress_every = suppressWarnings(as.integer(progress_every))
  if (length(progress_every) != 1L || is.na(progress_every) ||
      progress_every < 1L) {
    stop("`progress_every` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(progress_fun) && !is.function(progress_fun)) {
    stop("`progress_fun` must be a function or NULL.", call. = FALSE)
  }
  input = .tanimoto_normalize_input(
    compounds = compounds,
    compound_col = compound_col,
    cid_col = cid_col,
    inchikey_col = inchikey_col,
    smiles_col = smiles_col,
    compound_id_col = compound_id_col,
    group_cols = group_cols
  )
  compound_input = input$compounds
  membership = input$membership

  if (nrow(compound_input) < 2) {
    stop("At least two unique compounds are required for Tanimoto similarity.",
         call. = FALSE)
  }

  if (is.null(cache_dir)) cache_dir = .pubchem_default_cache_dir()
  if (isTRUE(refresh) && isTRUE(cache)) {
    unlink(file.path(cache_dir, "tanimoto"), recursive = TRUE, force = TRUE)
  }
  fetch = .pubchem_fetcher(cache = cache,
                           cache_dir = file.path(cache_dir, "tanimoto"),
                           throttle = throttle,
                           request_fun = request_fun,
                           service_busy_limit = service_busy_limit,
                           event_fun = progress_fun)

  resolution = .tanimoto_resolve_cids(compound_input, fetch,
                                      name_fallback = name_fallback,
                                      progress_fun = progress_fun,
                                      progress_every = progress_every)
  fingerprints = .tanimoto_fetch_fingerprints(
    resolution, fetch,
    progress_fun = progress_fun,
    progress_every = progress_every
  )
  canonical = .tanimoto_canonicalize_fingerprints(fingerprints, membership)
  fingerprints = canonical$fingerprints
  fingerprinted = canonical$compounds
  membership = canonical$membership
  if (nrow(fingerprinted) < 2) {
    stop(paste(
      "Fewer than two distinct, identity-consistent compounds resolved to",
      "PubChem Fingerprint2D values. Inspect `PubChemFingerprints` and",
      "`ExcludedCompoundIdentities` for identity mismatches or missing",
      "fingerprints."
    ), call. = FALSE)
  }

  bits = .tanimoto_decode_fingerprints(fingerprinted)
  fingerprinted = fingerprinted[match(rownames(bits), fingerprinted$compound_id),
                                , drop = FALSE]

  manifest = .tanimoto_empty_manifest()
  compound_pairs = .uaf_empty_table(.tanimoto_compound_pair_cols())
  if (isTRUE(return_compound_pairs)) {
    pair_count = choose(nrow(bits), 2)
    if (!is.null(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      path = file.path(out_dir, "compound_pair_tanimoto.csv.gz")
      written = .tanimoto_write_compound_pairs(
        bits, fingerprinted, path,
        block_size = pair_block_size,
        shard_rows = pair_shard_rows,
        progress_fun = progress_fun,
        progress_every = progress_every
      )
      manifest = rbind(manifest, .tanimoto_manifest_rows(
        "CompoundTanimoto", written,
        "Pairwise PubChem Fingerprint2D Tanimoto between unique compounds."
      ))
      compound_pairs = .uaf_empty_table(.tanimoto_compound_pair_cols())
    } else {
      if (pair_count > max_in_memory_pairs) {
        stop("Compound pair table would contain ", pair_count,
             " rows. Supply `out_dir` to stream the compressed CSV output, ",
             "or increase `max_in_memory_pairs` deliberately.",
             call. = FALSE)
      }
      compound_pairs = .tanimoto_compound_pair_table(bits, fingerprinted,
                                                     block_size = pair_block_size)
    }
  }

  group_summary = .uaf_empty_table(.tanimoto_group_summary_cols(group_cols))
  group_pairs = .uaf_empty_table(.tanimoto_group_pair_cols(group_cols))
  comparable_scope_summary = .uaf_empty_table(
    .comparable_pair_summary_cols("scope", thresholds)
  )
  comparable_group_summary = .uaf_empty_table(
    .comparable_pair_summary_cols("group", thresholds)
  )
  if (length(group_cols) > 0 && nrow(membership) > 0) {
    membership = membership[membership$compound_id %in% rownames(bits), ,
                            drop = FALSE]
    group_summary = .tanimoto_group_pair_summary(
      bits = bits,
      membership = membership,
      compounds = fingerprinted,
      group_cols = group_cols,
      thresholds = thresholds,
      top_n_pairs = top_n_pairs
    )
    comparable_scope_summary = .tanimoto_comparable_pair_summaries(
      bits, membership, fingerprinted, group_cols, thresholds, top_n_pairs,
      field = "comparison_scope", type = "scope"
    )
    comparable_group_summary = .tanimoto_comparable_pair_summaries(
      bits, membership, fingerprinted, group_cols, thresholds, top_n_pairs,
      field = "comparison_group", type = "group"
    )
    if (isTRUE(return_group_compound_pairs)) {
      group_pair_count = .tanimoto_group_pair_count(membership)
      if (!is.null(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        path = file.path(out_dir, "group_compound_pair_tanimoto.csv.gz")
        written = .tanimoto_write_group_compound_pairs(
          bits, membership, fingerprinted, path,
          block_size = pair_block_size,
          shard_rows = pair_shard_rows,
          progress_fun = progress_fun,
          progress_every = progress_every
        )
        manifest = rbind(manifest, .tanimoto_manifest_rows(
          "GroupCompoundTanimoto", written,
          "Cross-group compound-pair PubChem Fingerprint2D Tanimoto table."
        ))
      } else {
        if (group_pair_count > max_in_memory_pairs) {
          stop("Group compound pair table would contain ", group_pair_count,
               " rows. Supply `out_dir` to stream the compressed CSV output, ",
               "or increase `max_in_memory_pairs` deliberately.",
               call. = FALSE)
        }
        group_pairs = .tanimoto_group_compound_pair_table(
          bits = bits,
          membership = membership,
          compounds = fingerprinted,
          block_size = pair_block_size
        )
      }
    }
  }

  result = list(
    CompoundInput = compound_input,
    CompoundResolution = resolution,
    PubChemFingerprints = fingerprints,
    CompoundIdentityMap = canonical$identity_map,
    ExcludedCompoundIdentities = canonical$excluded,
    CanonicalCompounds = fingerprinted,
    GroupCompoundMembershipEvidence = canonical$membership_evidence,
    CompoundTanimoto = compound_pairs,
    GroupCompoundMembership = membership,
    GroupPairTanimotoSummary = group_summary,
    ComparableScopeTanimotoSummary = comparable_scope_summary,
    ComparableGroupTanimotoSummary = comparable_group_summary,
    GroupCompoundTanimoto = group_pairs,
    ExportManifest = manifest,
    Provenance = data.frame(
      Function = "chemicalTanimotoSimilarity",
      FingerprintSource = "PubChem Fingerprint2D",
      RetrievedAt = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      CacheEnabled = isTRUE(cache),
      CacheDir = cache_dir,
      stringsAsFactors = FALSE
    )
  )
  class(result) = c("uaf_tanimoto_similarity", "list")
  result
}

#' Plant-labeled PubChem fingerprint Tanimoto similarity
#'
#' @description
#' `plantChemicalTanimotoSimilarity()` converts a plant phytochemistry result
#' from `resolvePlantPhytochemistry()` or `runPlantPhytochemistryBatch()` into a
#' species-compound membership table, then calls
#' `chemicalTanimotoSimilarity()`. The returned tables preserve plant labels,
#' evidence tiers, plant-part context, and source identifiers so pairwise
#' chemistry outputs can be joined directly to downstream ecological,
#' phylogenetic, remediation, or trait analyses.
#'
#' @param plant_chemistry A plant phytochemistry result list containing
#' `PlantCompoundOccurrences` and `CompoundResolution`, or a data frame already
#' containing plant-compound rows.
#' @param level Plant grouping level. Defaults to `"species"`.
#' @param include_review_required Logical. If `TRUE`, carry
#' `CompoundIdentityReview` flags into the membership table when present.
#' @param ... Additional arguments passed to `chemicalTanimotoSimilarity()`,
#' including `out_dir`, `cache_dir`, `return_group_compound_pairs`, and
#' `request_fun`.
#'
#' @return A `"uaf_plant_tanimoto_similarity"` list containing the generic
#' Tanimoto tables plus plant-labeled aliases: `PlantCompoundMembership`,
#' `PlantPairTanimotoSummary`, and `PlantCompoundTanimoto`. Generic exact
#' comparable-scope and comparable-group summaries remain available as
#' `ComparableScopeTanimotoSummary` and `ComparableGroupTanimotoSummary`.
#'
#' @examples
#' \dontrun{
#' phyto = resolvePlantPhytochemistry(c("Salix nigra", "Zea mays"))
#' sim = plantChemicalTanimotoSimilarity(
#'   phyto,
#'   out_dir = "plant_tanimoto",
#'   return_group_compound_pairs = TRUE
#' )
#' sim$PlantPairTanimotoSummary
#' }
#'
#' @export
plantChemicalTanimotoSimilarity = function(plant_chemistry,
                                           level = "species",
                                           include_review_required = TRUE,
                                           ...) {
  table = .tanimoto_plant_membership_table(
    plant_chemistry,
    level = level,
    include_review_required = include_review_required
  )
  if (!level %in% names(table)) {
    stop("Plant chemistry input does not contain grouping level `", level,
         "`.", call. = FALSE)
  }
  result = chemicalTanimotoSimilarity(
    table,
    compound_col = "compound_name",
    cid_col = "CID",
    inchikey_col = "InChIKey",
    smiles_col = "SMILES",
    compound_id_col = "compound_id",
    group_cols = level,
    ...
  )

  result$PlantCompoundMembership = result$GroupCompoundMembership
  result$PlantPairTanimotoSummary = .tanimoto_label_plant_summary(
    result$GroupPairTanimotoSummary,
    level = level
  )
  result$PlantCompoundTanimoto = .tanimoto_label_plant_pairs(
    result$GroupCompoundTanimoto,
    level = level
  )
  result$Provenance$Function = "plantChemicalTanimotoSimilarity"
  class(result) = c("uaf_plant_tanimoto_similarity",
                    "uaf_tanimoto_similarity", "list")
  result
}

#' Prepare plant chemistry for structure-based Tanimoto analysis
#'
#' @description
#' Builds a conservative, offline handoff between species-first plant chemistry
#' discovery and PubChem Fingerprint2D Tanimoto analysis. Evidence rows remain
#' available for provenance, while the analysis membership is collapsed to one
#' row per plant and source-backed chemical identity. Repeated publications,
#' source records, and aliases are summarized rather than treated as distinct
#' compounds. Exact source-record identities are preferred when
#' `SourceCompoundIdentity` is available. Any inherited fingerprint columns are
#' deliberately removed. No fingerprints or pairwise similarities are computed.
#'
#' Valid full InChIKeys are preferred as identity keys. A positive PubChem CID
#' is used only when an InChIKey is unavailable. When an InChIKey is present,
#' any existing CID is retained only as reported audit metadata and is not
#' passed as trusted input to the next PubChem stage. This forces the downstream
#' workflow to resolve and verify the InChIKey with the collision-resistant
#' uafR cache.
#'
#' `DuplicateAudit` reports repeated evidence collapsed within one
#' plant-structure membership. `NameStructureAudit` separately reports
#' normalized compound labels associated with multiple exact structures. Those
#' structures remain distinct and require label-level review; they are never
#' collapsed merely because their names match.
#'
#' @param plant_chemistry A plant phytochemistry result containing
#' `PlantCompoundOccurrences` and `CompoundResolution`, or a data frame that
#' already combines plant labels with source-backed identity fields.
#' @param level Grouping column, normally `"species"`.
#' @param occurrence_status Allowed occurrence statuses when that column is
#' present. Defaults to direct or curated reported evidence.
#' @param analysis_ready If `TRUE`, require analysis-ready evidence when the
#' input contains an `analysis_ready` column.
#' @param min_confidence Minimum evidence confidence when a `confidence` column
#' is present.
#' @param include_review_required If `FALSE`, identity rows flagged for review
#' are excluded from the Tanimoto-ready membership and retained in the excluded
#' table.
#' @param out_dir Optional directory for the CSV/JSON server handoff bundle.
#' @param overwrite Logical. If `FALSE`, existing handoff files are protected.
#' @param strict Logical. If `TRUE`, stop when a validation check fails.
#'
#' @return A list with `PlantCompoundMembership`, `CompoundInput`,
#' `EvidenceRows`, `DuplicateAudit`, `NameStructureAudit`, `ExcludedRows`,
#' `ValidationSummary`, `Summary`, `ExportManifest`, and `Provenance`.
#'
#' @examples
#' \dontrun{
#' prepared = preparePlantTanimotoInput(
#'   phyto,
#'   out_dir = "plant_tanimoto_handoff"
#' )
#' prepared$PlantCompoundMembership
#' }
#'
#' @export
preparePlantTanimotoInput = function(
    plant_chemistry,
    level = "species",
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE,
    min_confidence = "medium",
    include_review_required = FALSE,
    out_dir = NULL,
    overwrite = FALSE,
    strict = TRUE) {
  if (length(level) != 1L || is.na(level) || !nzchar(level)) {
    stop("`level` must be one non-empty column name.", call. = FALSE)
  }
  level = as.character(level[[1]])
  min_confidence = tolower(.uaf_squish_text(min_confidence))
  if (length(min_confidence) != 1L || is.na(min_confidence) ||
      !min_confidence %in% c("unknown", "low", "medium", "high")) {
    stop("`min_confidence` must be one of unknown, low, medium, or high.",
         call. = FALSE)
  }
  rows = .tanimoto_preparation_rows(plant_chemistry, level)
  if (nrow(rows) < 1) {
    stop("No plant-compound rows were available for Tanimoto preparation.",
         call. = FALSE)
  }

  allowed_status = tolower(.uaf_non_empty(occurrence_status))
  status = tolower(.uaf_squish_text(rows$occurrence_status))
  status_known = !is.na(status) & status != ""
  status_ok = !rows$.occurrence_status_available |
    length(allowed_status) < 1 |
    (status_known & status %in% allowed_status)
  ready = tolower(.uaf_squish_text(rows$analysis_ready))
  ready_known = !is.na(ready) & ready != ""
  ready_ok = !isTRUE(analysis_ready) |
    !rows$.analysis_ready_available |
    (ready_known & ready %in% c("yes", "true", "1", "ready"))
  confidence_score = .plant_confidence_score(rows$confidence)
  minimum_score = .plant_confidence_score(min_confidence)[[1]]
  confidence_known = !is.na(.uaf_squish_text(rows$confidence)) &
    .uaf_squish_text(rows$confidence) != ""
  confidence_ok = !rows$.confidence_available |
    (confidence_known & confidence_score >= minimum_score)

  inchikey = toupper(.uaf_squish_text(rows$InChIKey))
  valid_inchikey = .tanimoto_valid_inchikey(inchikey)
  reported_cid = .uaf_squish_text(rows$reported_CID)
  normalized_cid = .tanimoto_normalize_cid(rows$CID)
  reported_valid_cid = !is.na(.tanimoto_normalize_cid(reported_cid))
  valid_cid = !is.na(normalized_cid)
  rows$compound_id = NA_character_
  rows$compound_id[valid_inchikey] =
    .tanimoto_sanitize_id(inchikey[valid_inchikey])
  cid_only = !valid_inchikey & valid_cid
  rows$compound_id[cid_only] = paste0("cid_", normalized_cid[cid_only])
  rows$InChIKey = ifelse(valid_inchikey, inchikey, NA_character_)
  rows$CID = ifelse(cid_only, normalized_cid, NA_character_)
  rows$reported_CID_usage = ifelse(
    valid_inchikey & reported_valid_cid,
    "audit_only_inchikey_will_be_resolved_and_verified",
    ifelse(cid_only, "direct_cid_requires_pubchem_verification",
           "no_usable_reported_cid")
  )
  rows$structure_identity_source = ifelse(
    valid_inchikey & rows$source_identity_exact %in% TRUE,
    "exact_source_record_inchikey",
    ifelse(valid_inchikey, "source_inchikey",
    ifelse(cid_only, "source_cid_without_inchikey", "no_usable_identity")
    )
  )
  structure_ok = !is.na(rows$compound_id) & rows$compound_id != ""
  group_ok = !is.na(rows[[level]]) & rows[[level]] != ""
  review_ok = isTRUE(include_review_required) |
    !(rows$identity_review_required %in% TRUE)
  rows$evidence_policy_pass = status_ok & ready_ok & confidence_ok
  rows$tanimoto_input_eligible = rows$evidence_policy_pass & group_ok &
    structure_ok & review_ok
  rows$preparation_exclusion_reason = vapply(seq_len(nrow(rows)), function(i) {
    reasons = character()
    if (!status_ok[[i]]) reasons = c(reasons, "occurrence_status_not_allowed")
    if (!ready_ok[[i]]) reasons = c(reasons, "not_analysis_ready")
    if (!confidence_ok[[i]]) reasons = c(reasons, "below_minimum_confidence")
    if (!group_ok[[i]]) reasons = c(reasons, "missing_group_value")
    if (!structure_ok[[i]]) reasons = c(reasons,
                                        "missing_valid_inchikey_or_cid")
    if (!review_ok[[i]]) reasons = c(reasons, "identity_review_required")
    if (length(reasons) < 1) NA_character_ else paste(reasons, collapse = "; ")
  }, character(1))

  eligible = rows[rows$tanimoto_input_eligible, , drop = FALSE]
  excluded = rows[!rows$tanimoto_input_eligible, , drop = FALSE]
  membership = .tanimoto_prepare_membership(eligible, level)
  compounds = .tanimoto_prepare_compounds(membership, level)
  duplicate_audit = .tanimoto_preparation_duplicates(eligible, level)
  name_structure_audit = .tanimoto_name_structure_audit(eligible, level)
  validation = .tanimoto_preparation_validation(
    rows, membership, compounds, excluded, level,
    include_review_required = include_review_required
  )
  summary = data.frame(
    input_evidence_rows = nrow(rows),
    eligible_evidence_rows = nrow(eligible),
    excluded_evidence_rows = nrow(excluded),
    plant_count = length(unique(.uaf_non_empty(rows[[level]]))),
    plants_with_tanimoto_ready_compounds =
      length(unique(.uaf_non_empty(membership[[level]]))),
    plant_compound_membership_count = nrow(membership),
    unique_structure_count = nrow(compounds),
    collapsed_repeated_evidence_rows = nrow(eligible) - nrow(membership),
    name_structure_ambiguity_count = nrow(name_structure_audit),
    identity_review_exclusion_count =
      sum(grepl("identity_review_required",
                excluded$preparation_exclusion_reason, fixed = TRUE),
          na.rm = TRUE),
    missing_structure_exclusion_count =
      sum(grepl("missing_valid_inchikey_or_cid",
                excluded$preparation_exclusion_reason, fixed = TRUE),
          na.rm = TRUE),
    exact_source_record_identity_count =
      sum(rows$source_identity_exact %in% TRUE),
    ignored_input_fingerprint_column_count =
      length(attr(rows, "ignored_input_fingerprint_columns", exact = TRUE)),
    validation_status = ifelse(any(validation$status == "fail"),
                               "fail", "pass"),
    stringsAsFactors = FALSE
  )
  provenance = data.frame(
    workflow = "preparePlantTanimotoInput",
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    package_version = tryCatch(as.character(utils::packageVersion("uafR")),
                               error = function(e) "development"),
    grouping_level = level,
    allowed_occurrence_status = paste(allowed_status, collapse = "; "),
    min_confidence = min_confidence,
    analysis_ready_required = isTRUE(analysis_ready),
    review_required_included = isTRUE(include_review_required),
    ignored_input_fingerprint_columns = paste(
      attr(rows, "ignored_input_fingerprint_columns", exact = TRUE),
      collapse = "; "
    ),
    identity_rule = paste(
      "Full source InChIKey preferred; positive source CID used only when",
      "InChIKey is unavailable. Existing fingerprints are not reused."
    ),
    stringsAsFactors = FALSE
  )

  result = list(
    PlantCompoundMembership = membership,
    CompoundInput = compounds,
    EvidenceRows = rows,
    DuplicateAudit = duplicate_audit,
    NameStructureAudit = name_structure_audit,
    ExcludedRows = excluded,
    ValidationSummary = validation,
    Summary = summary,
    ExportManifest = data.frame(),
    Provenance = provenance
  )
  class(result) = c("uaf_plant_tanimoto_input", "list")
  if (isTRUE(strict) && any(validation$status == "fail")) {
    failed = validation$check[validation$status == "fail"]
    stop("Tanimoto input validation failed: ", paste(failed, collapse = "; "),
         call. = FALSE)
  }
  if (!is.null(out_dir)) {
    result$ExportManifest = .tanimoto_write_preparation_bundle(
      result, out_dir, overwrite
    )
  }
  result
}

.tanimoto_preparation_rows = function(x, level) {
  if (is.data.frame(x)) {
    rows = as.data.frame(x, stringsAsFactors = FALSE)
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    occurrences = as.data.frame(x$PlantCompoundOccurrences,
                                stringsAsFactors = FALSE)
    resolution = if (is.data.frame(x$CompoundResolution)) {
      as.data.frame(x$CompoundResolution, stringsAsFactors = FALSE)
    } else {
      data.frame()
    }
    source_identity = .tanimoto_source_identity_table(x, resolution)
    rows = .tanimoto_apply_source_record_identity(occurrences,
                                                   source_identity)
    for (col in c("compound_name_clean", "compound_name", "resolved", "CID",
                  "InChIKey", "SMILES", "MolecularFormula",
                  "resolution_source", "identity_review_required")) {
      if (!col %in% names(resolution)) {
        resolution[[col]] = rep(NA, nrow(resolution))
      }
    }
    resolution$compound_name_clean =
      .uaf_squish_text(resolution$compound_name_clean)
    resolution$.structure_score =
      as.integer(resolution$resolved %in% TRUE) * 10L +
      as.integer(!is.na(.uaf_squish_text(resolution$InChIKey))) * 4L +
      as.integer(!is.na(.uaf_squish_text(resolution$CID))) * 2L +
      as.integer(!is.na(.uaf_squish_text(resolution$SMILES)))
    resolution = resolution[order(resolution$compound_name_clean,
                                  -resolution$.structure_score), , drop = FALSE]
    resolution = resolution[!duplicated(resolution$compound_name_clean), ,
                            drop = FALSE]
    idx = match(.uaf_squish_text(rows$compound_name_clean),
                resolution$compound_name_clean)
    resolution_matched = !is.na(idx) & resolution$resolved[idx] %in% TRUE
    for (col in c("CID", "InChIKey", "SMILES", "MolecularFormula",
                  "resolution_source", "resolved",
                  "identity_review_required")) {
      existing = if (col %in% names(rows)) rows[[col]] else
        rep(NA, nrow(rows))
      supplied = resolution[[col]][idx]
      if (identical(col, "identity_review_required")) {
        existing = .tanimoto_logical_value(existing)
        supplied = .tanimoto_logical_value(supplied)
        rows[[col]] = existing | (supplied & !rows$source_identity_exact)
      } else if (identical(col, "resolved")) {
        existing = existing %in% TRUE
        supplied = supplied %in% TRUE
        rows[[col]] = existing | supplied
      } else {
        rows[[col]] = .tanimoto_first_non_empty_vector(existing, supplied)
      }
    }
    rows$resolution_identity_matched = resolution_matched
    untrusted_occurrence_identity = !rows$source_identity_exact &
      !rows$resolution_identity_matched
    for (col in c("CID", "InChIKey", "SMILES", "MolecularFormula")) {
      rows[[col]][untrusted_occurrence_identity] = NA
    }
    rows = .tanimoto_apply_identity_review(rows, x$CompoundIdentityReview)
    rows = .tanimoto_apply_comparability(rows, x$ChemistryComparability)
  } else {
    stop(paste(
      "`plant_chemistry` must be a plant phytochemistry result containing",
      "PlantCompoundOccurrences, or a plant-compound data frame with identity",
      "fields."
    ), call. = FALSE)
  }
  if (!level %in% names(rows)) {
    stop("Plant chemistry input does not contain grouping level `", level,
         "`.", call. = FALSE)
  }
  rows$.occurrence_status_available = "occurrence_status" %in% names(rows)
  rows$.analysis_ready_available = "analysis_ready" %in% names(rows)
  rows$.confidence_available = "confidence" %in% names(rows)
  if (!"source_identity_exact" %in% names(rows)) {
    rows$source_identity_exact = FALSE
  }
  if (!"source_identity_match_status" %in% names(rows)) {
    rows$source_identity_match_status = "not_supplied"
  }
  if (!"resolution_identity_matched" %in% names(rows)) {
    rows$resolution_identity_matched = TRUE
  }

  aliases = list(
    compound_name = c("compound_name", "name", "chemical"),
    compound_name_clean = c("compound_name_clean", "query", "query_clean"),
    CID = c("source_CID", "source_cid", "CID", "cid", "pubchem_cid"),
    InChIKey = c("InChIKey_source", "source_InChIKey", "source_inchikey",
                  "InChIKey", "inchikey", "pubchem_inchikey"),
    SMILES = c("SMILES_source", "source_SMILES", "source_smiles", "SMILES",
               "smiles", "IsomericSMILES", "CanonicalSMILES"),
    MolecularFormula = c("MolecularFormula_source", "source_MolecularFormula",
                         "source_molecular_formula", "MolecularFormula",
                         "molecular_formula"),
    occurrence_status = c("occurrence_status"),
    analysis_ready = c("analysis_ready"),
    confidence = c("confidence"),
    identity_review_required = c("identity_review_required",
                                 "review_required")
  )
  for (target in names(aliases)) {
    if (!target %in% names(rows)) {
      hit = aliases[[target]][aliases[[target]] %in% names(rows)]
      rows[[target]] = if (length(hit) > 0) rows[[hit[[1]]]] else NA
    }
  }

  source_preference = list(
    CID = c("source_CID", "source_cid"),
    InChIKey = c("InChIKey_source", "source_InChIKey", "source_inchikey"),
    SMILES = c("SMILES_source", "source_SMILES", "source_smiles"),
    MolecularFormula = c("MolecularFormula_source", "source_MolecularFormula",
                         "source_molecular_formula")
  )
  for (target in names(source_preference)) {
    hit = source_preference[[target]][source_preference[[target]] %in%
                                        names(rows)]
    if (length(hit) > 0) {
      rows[[target]] = .tanimoto_first_non_empty_vector(
        rows[[hit[[1]]]], rows[[target]]
      )
    }
  }

  rows[[level]] = .uaf_squish_text(rows[[level]])
  rows$compound_name = .uaf_squish_text(rows$compound_name)
  rows$compound_name_clean = .uaf_squish_text(rows$compound_name_clean)
  missing_clean = is.na(rows$compound_name_clean)
  rows$compound_name_clean[missing_clean] =
    .plant_clean_compound(rows$compound_name[missing_clean])
  if (!"input_reported_CID" %in% names(rows)) {
    rows$input_reported_CID = .uaf_squish_text(rows$CID)
  }
  rows$reported_CID = .tanimoto_first_non_empty_vector(
    rows$input_reported_CID, rows$CID
  )
  rows$identity_review_required =
    .tanimoto_logical_value(rows$identity_review_required)
  rows$source_row_number = seq_len(nrow(rows))
  fingerprint_cols = grep("fingerprint", names(rows), ignore.case = TRUE,
                          value = TRUE)
  fingerprint_cols = setdiff(fingerprint_cols,
                             "ignored_input_fingerprint_columns")
  if (length(fingerprint_cols) > 0) rows[fingerprint_cols] = NULL
  attr(rows, "ignored_input_fingerprint_columns") = fingerprint_cols
  rows
}

.tanimoto_apply_comparability = function(rows, comparability) {
  if (!is.data.frame(comparability) || nrow(comparability) < 1 ||
      nrow(rows) < 1) {
    return(rows)
  }
  comparability = as.data.frame(comparability, stringsAsFactors = FALSE)
  key_cols = c("species", "compound_name_clean", "source_database",
               "source_record_id")
  for (col in key_cols) {
    if (!col %in% names(rows)) rows[[col]] = NA_character_
    if (!col %in% names(comparability)) comparability[[col]] = NA_character_
  }
  value_cols = c(
    "metabolism_domain", "biosynthetic_family", "chemical_behavior",
    "comparison_scope", "comparison_group", "comparison_subgroup",
    "comparability_confidence", "comparability_basis",
    "classification_source", "classification_source_table",
    "classification_source_field", "classification_source_value",
    "classification_source_confidence", "comparable_for_matrix",
    "comparison_caveat"
  )
  for (col in value_cols) {
    if (!col %in% names(comparability)) comparability[[col]] = NA_character_
  }
  make_key = function(x, cols) {
    do.call(paste, c(lapply(cols, function(col) {
      tolower(.uaf_squish_text(x[[col]]))
    }), sep = "\r"))
  }
  exact_key = make_key(comparability, key_cols)
  row_key = make_key(rows, key_cols)
  idx = match(row_key, exact_key)

  # Older occurrence tables may not have source-record identifiers. A
  # species/compound fallback is accepted only when all matching rows agree on
  # the complete comparability classification.
  unmatched = which(is.na(idx))
  if (length(unmatched) > 0) {
    fallback_cols = c("species", "compound_name_clean")
    comp_fallback = make_key(comparability, fallback_cols)
    row_fallback = make_key(rows, fallback_cols)
    groups = split(seq_len(nrow(comparability)), comp_fallback)
    for (i in unmatched) {
      candidates = groups[[row_fallback[[i]]]]
      if (length(candidates) < 1) next
      signatures = unique(make_key(comparability[candidates, , drop = FALSE],
                                   value_cols))
      if (length(signatures) == 1L) idx[[i]] = candidates[[1L]]
    }
  }
  for (col in value_cols) {
    supplied = comparability[[col]][idx]
    if (!col %in% names(rows)) {
      rows[[col]] = supplied
    } else {
      rows[[col]] = .tanimoto_first_non_empty_vector(rows[[col]], supplied)
    }
  }
  rows
}

.tanimoto_source_identity_table = function(x, resolution) {
  candidates = list(
    x$SourceCompoundIdentity,
    if (is.list(x$CategorateResult)) {
      x$CategorateResult$SourceCompoundIdentity
    } else NULL,
    attr(x$CompoundResolution, "SourceCompoundIdentity", exact = TRUE),
    if (is.list(attr(x$CompoundResolution, "CategorateResult", exact = TRUE))) {
      attr(x$CompoundResolution,
           "CategorateResult", exact = TRUE)$SourceCompoundIdentity
    } else NULL,
    attr(resolution, "SourceCompoundIdentity", exact = TRUE),
    if (is.list(attr(resolution, "CategorateResult", exact = TRUE))) {
      attr(resolution,
           "CategorateResult", exact = TRUE)$SourceCompoundIdentity
    } else NULL
  )
  candidates = candidates[vapply(candidates, is.data.frame, logical(1))]
  candidates = candidates[vapply(candidates, nrow, integer(1)) > 0L]
  if (length(candidates) < 1) return(data.frame())
  out = candidates[[1]]
  required = c(
    "compound_name", "compound_name_clean", "source_database",
    "source_record_id", "source_compound_id", "source_compound_id_type",
    "CID", "InChIKey", "SMILES", "MolecularFormula", "evidence_url",
    "evidence_text", "identity_status", "identity_note"
  )
  for (col in setdiff(required, names(out))) out[[col]] = NA_character_
  out[, required, drop = FALSE]
}

.tanimoto_apply_source_record_identity = function(rows, source_identity) {
  rows = as.data.frame(rows, stringsAsFactors = FALSE)
  rows$input_reported_CID = if ("CID" %in% names(rows)) {
    .uaf_squish_text(rows$CID)
  } else if ("pubchem_cid" %in% names(rows)) {
    .uaf_squish_text(rows$pubchem_cid)
  } else {
    rep(NA_character_, nrow(rows))
  }
  identity_cols = c("CID", "InChIKey", "SMILES", "MolecularFormula",
                    "resolution_source", "resolved",
                    "identity_review_required")
  for (col in identity_cols) {
    if (!col %in% names(rows)) rows[[col]] = rep(NA, nrow(rows))
  }
  rows$source_identity_exact = FALSE
  rows$source_identity_match_status = "no_source_identity_table"
  rows$source_identity_row_count = 0L
  rows$source_identity_record_ids = NA_character_
  rows$source_identity_statuses = NA_character_
  rows$source_identity_reported_CIDs = NA_character_
  if (!is.data.frame(source_identity) || nrow(source_identity) < 1 ||
      nrow(rows) < 1) {
    return(rows)
  }

  for (col in c("source_database", "source_record_id", "compound_id")) {
    if (!col %in% names(rows)) rows[[col]] = NA_character_
  }
  source_identity = as.data.frame(source_identity, stringsAsFactors = FALSE)
  for (col in c("source_database", "source_record_id", "source_compound_id",
                "CID", "InChIKey", "SMILES", "MolecularFormula",
                "identity_status")) {
    if (!col %in% names(source_identity)) {
      source_identity[[col]] = rep(NA_character_, nrow(source_identity))
    }
  }

  source_db = tolower(.uaf_squish_text(source_identity$source_database))
  source_ids = lapply(seq_len(nrow(source_identity)), function(i) {
    unique(.uaf_non_empty(c(source_identity$source_record_id[[i]],
                            source_identity$source_compound_id[[i]])))
  })
  map_keys = unlist(lapply(seq_len(nrow(source_identity)), function(i) {
    ids = source_ids[[i]]
    if (length(ids) < 1) return(character())
    db = source_db[[i]]
    unique(c(paste0("*\r", ids),
             if (!is.na(db) && db != "") paste0(db, "\r", ids) else
               character()))
  }), use.names = FALSE)
  map_rows = unlist(lapply(seq_len(nrow(source_identity)), function(i) {
    ids = source_ids[[i]]
    if (length(ids) < 1) return(integer())
    db = source_db[[i]]
    rep(i, length(ids) * ifelse(!is.na(db) && db != "", 2L, 1L))
  }), use.names = FALSE)
  source_map = split(map_rows, map_keys)

  for (i in seq_len(nrow(rows))) {
    ids = unique(.uaf_non_empty(c(rows$source_record_id[[i]],
                                  rows$compound_id[[i]])))
    if (length(ids) < 1) {
      rows$source_identity_match_status[[i]] = "no_source_record_identifier"
      next
    }
    db = tolower(.uaf_squish_text(rows$source_database[[i]]))
    keys = paste0("*\r", ids)
    if (!is.na(db) && db != "") keys = c(paste0(db, "\r", ids), keys)
    idx = unique(unlist(source_map[keys], use.names = FALSE))
    if (length(idx) < 1) {
      rows$source_identity_match_status[[i]] = "source_record_not_matched"
      next
    }
    part = source_identity[idx, , drop = FALSE]
    part_db = tolower(.uaf_squish_text(part$source_database))
    if (!is.na(db) && db != "") {
      same_db = is.na(part_db) | part_db == "" | part_db == db
      part = part[same_db, , drop = FALSE]
    }
    if (nrow(part) < 1) {
      rows$source_identity_match_status[[i]] = "source_database_mismatch"
      next
    }

    keys_found = toupper(.uaf_squish_text(part$InChIKey))
    keys_found = unique(keys_found[.tanimoto_valid_inchikey(keys_found)])
    cids_found = .tanimoto_normalize_cid(part$CID)
    cids_found = unique(.uaf_non_empty(cids_found))
    rows$source_identity_row_count[[i]] = nrow(part)
    rows$source_identity_record_ids[[i]] = .tanimoto_collapse_values(c(
      part$source_record_id, part$source_compound_id
    ))
    rows$source_identity_statuses[[i]] =
      .tanimoto_collapse_values(part$identity_status)
    rows$source_identity_reported_CIDs[[i]] =
      .tanimoto_collapse_values(cids_found)

    exact_key = length(keys_found) == 1L
    exact_cid = length(keys_found) < 1L && length(cids_found) == 1L
    if (!exact_key && !exact_cid) {
      rows$source_identity_match_status[[i]] = if (
        length(keys_found) > 1L || length(cids_found) > 1L
      ) "source_record_identity_conflict" else "source_record_structure_missing"
      rows$identity_review_required[[i]] = TRUE
      next
    }

    selected = if (exact_key) {
      part[toupper(.uaf_squish_text(part$InChIKey)) == keys_found[[1]], ,
           drop = FALSE]
    } else {
      part[.tanimoto_normalize_cid(part$CID) == cids_found[[1]], , drop = FALSE]
    }
    selected = selected[order(.uaf_squish_text(selected$source_record_id),
                              .uaf_squish_text(selected$source_compound_id),
                              na.last = TRUE), , drop = FALSE]
    selected = selected[1, , drop = FALSE]
    rows$InChIKey[[i]] = if (exact_key) keys_found[[1]] else NA_character_
    rows$CID[[i]] = if (exact_cid) cids_found[[1]] else
      .uaf_first_non_empty_text(selected$CID)
    rows$SMILES[[i]] = .uaf_first_non_empty_text(selected$SMILES)
    rows$MolecularFormula[[i]] =
      .uaf_first_non_empty_text(selected$MolecularFormula)
    rows$resolution_source[[i]] = paste0(
      toupper(.uaf_first_non_empty_text(selected$source_database, "source")),
      "_exact_source_record_identity"
    )
    rows$resolved[[i]] = TRUE
    rows$source_identity_exact[[i]] = TRUE
    rows$source_identity_match_status[[i]] = if (exact_key) {
      "exact_source_record_inchikey"
    } else {
      "exact_source_record_cid"
    }
  }
  rows
}

.tanimoto_apply_identity_review = function(rows, review) {
  if (!is.data.frame(review) || nrow(review) < 1 || nrow(rows) < 1) {
    return(rows)
  }
  for (col in c("compound_name_clean", "identity_issue_type",
                "review_required", "review_decision")) {
    if (!col %in% names(review)) review[[col]] = NA_character_
  }
  review$compound_name_clean = .uaf_squish_text(review$compound_name_clean)
  review = review[!is.na(review$compound_name_clean) &
                    review$compound_name_clean != "", , drop = FALSE]
  if (nrow(review) < 1) return(rows)
  groups = split(seq_len(nrow(review)), review$compound_name_clean)
  for (key in names(groups)) {
    part = review[groups[[key]], , drop = FALSE]
    row_idx = which(.uaf_squish_text(rows$compound_name_clean) == key)
    if (length(row_idx) < 1) next
    decisions = tolower(.uaf_squish_text(part$review_decision))
    accepted = decisions %in% c(
      "accept", "accepted", "accept_resolved", "accept_resolved_identity",
      "accept_identity", "keep_resolved", "update", "update_identity",
      "replace", "replace_identity", "manual_identity", "curate_identity",
      "promote_identity"
    )
    explicitly_required = .tanimoto_logical_value(part$review_required)
    issue_types = unique(.uaf_non_empty(part$identity_issue_type))
    requires_review = if (all(accepted)) FALSE else
      any(explicitly_required) || length(issue_types) > 0L || nrow(part) > 0L
    only_name_level_source_ambiguity = length(issue_types) > 0L &&
      all(issue_types == "resolved_source_structure_ambiguous")
    per_row_required = rep(requires_review, length(row_idx))
    per_row_required[
      rows$source_identity_exact[row_idx] & only_name_level_source_ambiguity
    ] = FALSE
    rows$identity_review_required[row_idx] =
      .tanimoto_logical_value(rows$identity_review_required[row_idx]) |
      per_row_required
    if (!"identity_review_issue_types" %in% names(rows)) {
      rows$identity_review_issue_types = NA_character_
    }
    rows$identity_review_issue_types[row_idx] =
      .tanimoto_collapse_values(issue_types)
  }
  rows
}

.tanimoto_logical_value = function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  if (is.numeric(x)) return(!is.na(x) & x != 0)
  tolower(.uaf_squish_text(x)) %in%
    c("true", "yes", "y", "1", "review_required", "required")
}

.tanimoto_valid_inchikey = function(x) {
  x = toupper(.uaf_squish_text(x))
  !is.na(x) & grepl("^[A-Z]{14}-[A-Z]{10}-[A-Z]$", x)
}

.tanimoto_normalize_cid = function(x) {
  x = .uaf_squish_text(x)
  valid = !is.na(x) & grepl("^[0-9]+$", x)
  numeric = suppressWarnings(as.numeric(x))
  valid = valid & is.finite(numeric) & numeric > 0 & numeric == floor(numeric)
  out = rep(NA_character_, length(x))
  out[valid] = format(numeric[valid], scientific = FALSE, trim = TRUE)
  out
}

.tanimoto_first_non_empty_vector = function(primary, fallback) {
  primary = .uaf_squish_text(primary)
  fallback = .uaf_squish_text(fallback)
  missing = is.na(primary) | primary == ""
  primary[missing] = fallback[missing]
  primary
}

.tanimoto_prepare_membership = function(rows, level) {
  if (nrow(rows) < 1) return(.tanimoto_empty_prepared_membership(level))
  key = paste(rows[[level]], rows$compound_id, sep = "\r")
  groups = split(seq_len(nrow(rows)), key)
  out = lapply(groups, function(idx) {
    part = rows[idx, , drop = FALSE]
    names_observed = .uaf_non_empty(part$compound_name)
    if (length(names_observed) > 0) {
      name_counts = sort(table(names_observed), decreasing = TRUE)
      best_names = names(name_counts)[name_counts == max(name_counts)]
      compound_name = sort(best_names)[[1]]
    } else {
      compound_name = NA_character_
    }
    inchikeys = unique(.uaf_non_empty(part$InChIKey))
    cids = unique(.uaf_non_empty(part$CID))
    smiles = unique(.uaf_non_empty(part$SMILES))
    formulas = unique(.uaf_non_empty(part$MolecularFormula))
    comparison_scopes = unique(.uaf_non_empty(part$comparison_scope))
    comparison_groups = unique(.uaf_non_empty(part$comparison_group))
    comparison_subgroups = unique(.uaf_non_empty(part$comparison_subgroup))
    comparable_values = tolower(.uaf_squish_text(
      part$comparable_for_matrix
    ))
    comparable = length(comparison_scopes) == 1L &&
      length(comparison_groups) == 1L &&
      !tolower(comparison_scopes) %in% c("unknown", "mixed") &&
      !tolower(comparison_groups) %in% c("unknown", "mixed") &&
      length(comparable_values) > 0L &&
      all(comparable_values %in% c("yes", "true", "1"))
    data.frame(
      group_value = part[[level]][[1]],
      compound_id = part$compound_id[[1]],
      compound_name = compound_name,
      compound_aliases = .tanimoto_collapse_values(part$compound_name),
      InChIKey = if (length(inchikeys) == 1L) inchikeys else NA_character_,
      CID = if (length(inchikeys) < 1L && length(cids) == 1L) cids else
        NA_character_,
      reported_CIDs = .tanimoto_collapse_values(part$reported_CID),
      SMILES = if (length(smiles) == 1L) smiles else NA_character_,
      MolecularFormula = if (length(formulas) == 1L) formulas else
        NA_character_,
      structure_identity_source = part$structure_identity_source[[1]],
      source_identity_exact = any(part$source_identity_exact %in% TRUE),
      source_identity_match_statuses =
        .tanimoto_collapse_values(part$source_identity_match_status),
      evidence_row_count = nrow(part),
      source_record_count = if ("source_record_id" %in% names(part)) {
        length(unique(.uaf_non_empty(part$source_record_id)))
      } else NA_integer_,
      source_record_ids = if ("source_record_id" %in% names(part)) {
        .tanimoto_collapse_values(part$source_record_id)
      } else NA_character_,
      source_database_count = if ("source_database" %in% names(part)) {
        length(unique(.uaf_non_empty(part$source_database)))
      } else NA_integer_,
      source_databases = if ("source_database" %in% names(part)) {
        .tanimoto_collapse_values(part$source_database)
      } else NA_character_,
      evidence_tiers = if ("evidence_tier" %in% names(part)) {
        .tanimoto_collapse_values(part$evidence_tier)
      } else NA_character_,
      confidence_levels = .tanimoto_collapse_values(part$confidence),
      occurrence_statuses =
        .tanimoto_collapse_values(part$occurrence_status),
      plant_part_groups = if ("plant_part_group" %in% names(part)) {
        .tanimoto_collapse_values(part$plant_part_group)
      } else NA_character_,
      tissue_groups = if ("tissue_group" %in% names(part)) {
        .tanimoto_collapse_values(part$tissue_group)
      } else NA_character_,
      method_groups = if ("method_group" %in% names(part)) {
        .tanimoto_collapse_values(part$method_group)
      } else NA_character_,
      metabolism_domain = .tanimoto_single_classification(
        part$metabolism_domain
      ),
      biosynthetic_family = .tanimoto_single_classification(
        part$biosynthetic_family
      ),
      chemical_behavior = .tanimoto_single_classification(
        part$chemical_behavior
      ),
      comparison_scope = if (length(comparison_scopes) == 1L) {
        comparison_scopes[[1L]]
      } else if (length(comparison_scopes) > 1L) "mixed" else "unknown",
      comparison_group = if (length(comparison_groups) == 1L) {
        comparison_groups[[1L]]
      } else if (length(comparison_groups) > 1L) "mixed" else "unknown",
      comparison_subgroup = if (length(comparison_subgroups) == 1L) {
        comparison_subgroups[[1L]]
      } else if (length(comparison_subgroups) > 1L) "mixed" else "unknown",
      comparability_confidence = .tanimoto_single_classification(
        part$comparability_confidence
      ),
      comparability_basis = .tanimoto_collapse_values(
        part$comparability_basis
      ),
      classification_source = .tanimoto_single_classification(
        part$classification_source
      ),
      classification_sources = .tanimoto_collapse_values(
        part$classification_source
      ),
      comparable_for_matrix = ifelse(comparable, "Yes", "No"),
      comparison_caveats = .tanimoto_collapse_values(part$comparison_caveat),
      identity_review_required = any(part$identity_review_required %in% TRUE),
      tanimoto_input_status = if (
        any(part$identity_review_required %in% TRUE)
      ) "ready_with_explicit_review_override" else
        "ready_for_pubchem_fingerprint_resolution",
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, out)
  names(out)[names(out) == "group_value"] = level
  row.names(out) = NULL
  out[order(out[[level]], out$compound_id), , drop = FALSE]
}

.tanimoto_single_classification = function(x) {
  values = unique(.uaf_non_empty(x))
  if (length(values) < 1) return("unknown")
  if (length(values) == 1L) return(values[[1L]])
  "mixed"
}

.tanimoto_empty_prepared_membership = function(level) {
  out = data.frame(
    group_value = character(), compound_id = character(),
    compound_name = character(), compound_aliases = character(),
    InChIKey = character(), CID = character(), reported_CIDs = character(),
    SMILES = character(), MolecularFormula = character(),
    structure_identity_source = character(), source_identity_exact = logical(),
    source_identity_match_statuses = character(),
    evidence_row_count = integer(),
    source_record_count = integer(), source_record_ids = character(),
    source_database_count = integer(), source_databases = character(),
    evidence_tiers = character(), confidence_levels = character(),
    occurrence_statuses = character(), plant_part_groups = character(),
    tissue_groups = character(), method_groups = character(),
    metabolism_domain = character(), biosynthetic_family = character(),
    chemical_behavior = character(), comparison_scope = character(),
    comparison_group = character(), comparison_subgroup = character(),
    comparability_confidence = character(), comparability_basis = character(),
    classification_source = character(), classification_sources = character(),
    comparable_for_matrix = character(), comparison_caveats = character(),
    identity_review_required = logical(), tanimoto_input_status = character(),
    stringsAsFactors = FALSE
  )
  names(out)[names(out) == "group_value"] = level
  out
}

.tanimoto_prepare_compounds = function(membership, level) {
  if (nrow(membership) < 1) {
    return(data.frame(
      compound_id = character(), compound_name = character(),
      compound_aliases = character(), InChIKey = character(),
      CID = character(), reported_CIDs = character(), SMILES = character(),
      MolecularFormula = character(), structure_identity_source = character(),
      source_identity_exact = logical(),
      source_identity_match_statuses = character(), plant_count = integer(),
      evidence_row_count = integer(),
      stringsAsFactors = FALSE
    ))
  }
  groups = split(seq_len(nrow(membership)), membership$compound_id)
  out = lapply(groups, function(idx) {
    part = membership[idx, , drop = FALSE]
    names_observed = sort(unique(.uaf_non_empty(part$compound_name)))
    data.frame(
      compound_id = part$compound_id[[1]],
      compound_name = if (length(names_observed) > 0) {
        names_observed[[1]]
      } else {
        NA_character_
      },
      compound_aliases = .tanimoto_collapse_values(c(
        part$compound_name, unlist(strsplit(
          .uaf_non_empty(part$compound_aliases), "; ", fixed = TRUE
        ), use.names = FALSE)
      )),
      InChIKey = .uaf_first_non_empty_text(part$InChIKey),
      CID = .uaf_first_non_empty_text(part$CID),
      reported_CIDs = .tanimoto_collapse_values(part$reported_CIDs),
      SMILES = .uaf_first_non_empty_text(part$SMILES),
      MolecularFormula = .uaf_first_non_empty_text(part$MolecularFormula),
      structure_identity_source =
        .uaf_first_non_empty_text(part$structure_identity_source),
      source_identity_exact = any(part$source_identity_exact %in% TRUE),
      source_identity_match_statuses =
        .tanimoto_collapse_values(part$source_identity_match_statuses),
      plant_count = length(unique(.uaf_non_empty(part[[level]]))),
      evidence_row_count = sum(part$evidence_row_count, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, out)
  row.names(out) = NULL
  out[order(out$compound_id), , drop = FALSE]
}

.tanimoto_preparation_duplicates = function(rows, level) {
  if (nrow(rows) < 1) {
    out = data.frame(
      group_value = character(), compound_id = character(),
      evidence_row_count = integer(), compound_alias_count = integer(),
      compound_aliases = character(), source_record_count = integer(),
      duplicate_interpretation = character(), stringsAsFactors = FALSE
    )
    names(out)[names(out) == "group_value"] = level
    return(out)
  }
  key = paste(rows[[level]], rows$compound_id, sep = "\r")
  groups = split(seq_len(nrow(rows)), key)
  groups = groups[vapply(groups, length, integer(1)) > 1L]
  if (length(groups) < 1) {
    out = data.frame(
      group_value = character(), compound_id = character(),
      evidence_row_count = integer(), compound_alias_count = integer(),
      compound_aliases = character(), source_record_count = integer(),
      duplicate_interpretation = character(), stringsAsFactors = FALSE
    )
    names(out)[names(out) == "group_value"] = level
    return(out)
  }
  out = do.call(rbind, lapply(groups, function(idx) {
    part = rows[idx, , drop = FALSE]
    data.frame(
      group_value = part[[level]][[1]],
      compound_id = part$compound_id[[1]],
      evidence_row_count = nrow(part),
      compound_alias_count = length(unique(.uaf_non_empty(
        part$compound_name
      ))),
      compound_aliases = .tanimoto_collapse_values(part$compound_name),
      source_record_count = if ("source_record_id" %in% names(part)) {
        length(unique(.uaf_non_empty(part$source_record_id)))
      } else NA_integer_,
      duplicate_interpretation = paste(
        "Repeated evidence and/or aliases collapsed to one plant-structure",
        "membership; source rows remain in EvidenceRows."
      ),
      stringsAsFactors = FALSE
    )
  }))
  names(out)[names(out) == "group_value"] = level
  row.names(out) = NULL
  out[order(out[[level]], out$compound_id), , drop = FALSE]
}

.tanimoto_name_structure_audit = function(rows, level) {
  empty = data.frame(
    audit_scope = character(), group_value = character(),
    compound_name_clean = character(), compound_names = character(),
    structure_count = integer(), compound_ids = character(),
    group_count = integer(), groups = character(),
    evidence_row_count = integer(), source_record_count = integer(),
    source_record_ids = character(), source_databases = character(),
    interpretation = character(), stringsAsFactors = FALSE
  )
  names(empty)[names(empty) == "group_value"] = level
  if (nrow(rows) < 1) return(empty)

  clean_name = .uaf_squish_text(rows$compound_name_clean)
  missing_name = is.na(clean_name) | clean_name == ""
  clean_name[missing_name] = .plant_clean_compound(
    rows$compound_name[missing_name]
  )
  usable = !is.na(clean_name) & clean_name != "" &
    !is.na(rows$compound_id) & rows$compound_id != ""
  if (!any(usable)) return(empty)
  audit_rows = rows[usable, , drop = FALSE]
  audit_rows$.audit_compound_name_clean = clean_name[usable]

  summarize_groups = function(index_groups, scope) {
    index_groups = index_groups[vapply(index_groups, function(idx) {
      length(unique(.uaf_non_empty(audit_rows$compound_id[idx]))) > 1L
    }, logical(1))]
    if (length(index_groups) < 1) return(NULL)
    do.call(rbind, lapply(index_groups, function(idx) {
      part = audit_rows[idx, , drop = FALSE]
      group_values = sort(unique(.uaf_non_empty(part[[level]])))
      source_record_ids = if ("source_record_id" %in% names(part)) {
        .uaf_non_empty(part$source_record_id)
      } else character()
      source_databases = if ("source_database" %in% names(part)) {
        part$source_database
      } else character()
      data.frame(
        audit_scope = scope,
        group_value = if (identical(scope, "within_group")) {
          group_values[[1]]
        } else {
          NA_character_
        },
        compound_name_clean = part$.audit_compound_name_clean[[1]],
        compound_names = .tanimoto_collapse_values(part$compound_name),
        structure_count = length(unique(.uaf_non_empty(part$compound_id))),
        compound_ids = .tanimoto_collapse_values(part$compound_id),
        group_count = length(group_values),
        groups = paste(group_values, collapse = "; "),
        evidence_row_count = nrow(part),
        source_record_count = length(unique(source_record_ids)),
        source_record_ids = paste(sort(unique(source_record_ids)),
                                  collapse = "; "),
        source_databases = .tanimoto_collapse_values(source_databases),
        interpretation = if (identical(scope, "within_group")) {
          paste(
            "One normalized compound label is tied to multiple exact",
            "source-backed structures in the same plant. Structures remain",
            "separate and require label-level review before interpretation."
          )
        } else {
          paste(
            "One normalized compound label maps to multiple exact",
            "source-backed structures across the dataset. Structures remain",
            "separate and must not be collapsed by compound name."
          )
        },
        stringsAsFactors = FALSE
      )
    }))
  }

  global_groups = split(
    seq_len(nrow(audit_rows)),
    audit_rows$.audit_compound_name_clean,
    drop = TRUE
  )
  within_groups = split(
    seq_len(nrow(audit_rows)),
    paste(audit_rows[[level]], audit_rows$.audit_compound_name_clean,
          sep = "\r"),
    drop = TRUE
  )
  out = rbind(
    summarize_groups(global_groups, "global"),
    summarize_groups(within_groups, "within_group")
  )
  if (is.null(out) || nrow(out) < 1) return(empty)
  names(out)[names(out) == "group_value"] = level
  row.names(out) = NULL
  out[order(out$audit_scope, out$compound_name_clean, out[[level]],
            na.last = TRUE), , drop = FALSE]
}

.tanimoto_preparation_validation = function(rows, membership, compounds,
                                             excluded, level,
                                             include_review_required = FALSE) {
  membership_key = paste(membership[[level]], membership$compound_id,
                         sep = "\r")
  checks = data.frame(
    check = c(
      "input_rows_accounted_for", "membership_nonempty",
      "membership_group_nonmissing", "membership_key_unique",
      "membership_identity_nonmissing", "membership_identity_format_valid",
      "membership_review_flags_excluded", "compound_table_key_unique",
      "compound_table_matches_membership", "fingerprint_payload_absent",
      "at_least_two_structures"
    ),
    status = c(
      ifelse(nrow(rows) == sum(rows$tanimoto_input_eligible) + nrow(excluded),
             "pass", "fail"),
      ifelse(nrow(membership) > 0, "pass", "fail"),
      ifelse(all(!is.na(membership[[level]]) & membership[[level]] != ""),
             "pass", "fail"),
      ifelse(anyDuplicated(membership_key) == 0L, "pass", "fail"),
      ifelse(all(!is.na(membership$compound_id) &
                   (!is.na(membership$InChIKey) |
                      !is.na(membership$CID))), "pass", "fail"),
      ifelse(all(
        (.tanimoto_valid_inchikey(membership$InChIKey) |
           !is.na(.tanimoto_normalize_cid(membership$CID))) &
          !( !is.na(membership$InChIKey) & !is.na(membership$CID) )
      ), "pass", "fail"),
      ifelse(isTRUE(include_review_required) ||
               !any(membership$identity_review_required %in% TRUE),
             "pass", "fail"),
      ifelse(anyDuplicated(compounds$compound_id) == 0L, "pass", "fail"),
      ifelse(setequal(compounds$compound_id,
                      unique(membership$compound_id)), "pass", "fail"),
      ifelse(!any(grepl("fingerprint", c(names(membership), names(compounds)),
                        ignore.case = TRUE)), "pass", "fail"),
      ifelse(nrow(compounds) >= 2L, "pass", "fail")
    ),
    detail = c(
      paste(nrow(rows), "input rows;", sum(rows$tanimoto_input_eligible),
            "eligible;", nrow(excluded), "excluded"),
      paste(nrow(membership), "plant-structure memberships"),
      paste(sum(is.na(membership[[level]]) | membership[[level]] == ""),
            "memberships lack a grouping value"),
      paste(length(unique(membership_key)), "unique membership keys"),
      "Every prepared membership requires a valid InChIKey or source CID.",
      paste(
        "Each membership has exactly one valid full InChIKey or positive CID;",
        "CID is blank when an InChIKey is supplied."
      ),
      paste(sum(membership$identity_review_required %in% TRUE),
            "review-required prepared rows"),
      paste(nrow(compounds), "unique structure rows"),
      "CompoundInput and membership structure identifiers must agree.",
      "Prepared analysis tables contain no inherited fingerprint payload.",
      paste(nrow(compounds), "structures available")
    ),
    stringsAsFactors = FALSE
  )
  checks
}

.tanimoto_write_preparation_bundle = function(result, out_dir, overwrite) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tables = list(
    plant_compound_membership_tanimoto_ready =
      result$PlantCompoundMembership,
    compound_input_tanimoto_ready = result$CompoundInput,
    plant_compound_evidence_rows = result$EvidenceRows,
    plant_compound_duplicate_audit = result$DuplicateAudit,
    plant_compound_name_structure_audit = result$NameStructureAudit,
    plant_compound_excluded_from_tanimoto = result$ExcludedRows,
    tanimoto_input_validation = result$ValidationSummary,
    tanimoto_input_summary = result$Summary,
    tanimoto_input_provenance = result$Provenance
  )
  table_paths = file.path(out_dir, paste0(names(tables), ".csv"))
  manifest_path = file.path(out_dir, "tanimoto_input_export_manifest.csv")
  json_path = file.path(out_dir, "tanimoto_input_manifest.json")
  protected_paths = c(table_paths, manifest_path, json_path)
  existing = protected_paths[file.exists(protected_paths)]
  if (length(existing) > 0 && !isTRUE(overwrite)) {
    stop(
      "Tanimoto handoff files exist and `overwrite = FALSE`: ",
      paste(basename(existing), collapse = ", "),
      call. = FALSE
    )
  }
  manifest = lapply(names(tables), function(name) {
    path = file.path(out_dir, paste0(name, ".csv"))
    .plant_atomic_write_csv(tables[[name]], path)
    data.frame(
      artifact = name,
      file = basename(path),
      row_count = nrow(tables[[name]]),
      column_count = ncol(tables[[name]]),
      bytes = file.info(path)$size,
      md5 = unname(tools::md5sum(path)[[1]]),
      stringsAsFactors = FALSE
    )
  })
  manifest = do.call(rbind, manifest)
  .plant_atomic_write_csv(manifest, manifest_path)
  temp_json = tempfile(paste0(basename(json_path), "."),
                       tmpdir = dirname(json_path))
  on.exit(unlink(temp_json, force = TRUE), add = TRUE)
  jsonlite::write_json(
    list(
      workflow = "preparePlantTanimotoInput",
      schema_version = "1.0.0",
      package_version = result$Provenance$package_version[[1]],
      created_at = result$Provenance$created_at[[1]],
      summary = result$Summary,
      validation = result$ValidationSummary,
      files = manifest,
      csv_manifest = list(
        file = basename(manifest_path),
        bytes = unname(file.info(manifest_path)$size),
        md5 = unname(tools::md5sum(manifest_path)[[1]])
      ),
      next_step = paste(
        "Run plantChemicalTanimotoSimilarity() or",
        "chemicalTanimotoSimilarity(group_cols = 'species') on a server",
        "using plant_compound_membership_tanimoto_ready.csv."
      )
    ),
    temp_json, pretty = TRUE, auto_unbox = TRUE, na = "null",
    dataframe = "rows"
  )
  if (!file.rename(temp_json, json_path)) {
    stop("Could not atomically write: ", json_path, call. = FALSE)
  }
  manifest
}

.tanimoto_normalize_input = function(compounds, compound_col, cid_col,
                                     inchikey_col, smiles_col,
                                     compound_id_col, group_cols) {
  if (is.character(compounds)) {
    raw = data.frame(compound_name = compounds, stringsAsFactors = FALSE)
  } else if (is.data.frame(compounds)) {
    raw = compounds
  } else {
    stop("`compounds` must be a character vector or data frame.",
         call. = FALSE)
  }
  if (nrow(raw) < 1) stop("`compounds` contains no rows.", call. = FALSE)

  compound_col = .tanimoto_detect_col(raw, compound_col,
                                      c("compound_name", "Query", "Chemical",
                                        "name", "compound"))
  cid_col = .tanimoto_detect_col(raw, cid_col,
                                 c("CID", "pubchem_cid", "cid",
                                   "PubChemCID"))
  inchikey_col = .tanimoto_detect_col(raw, inchikey_col,
                                      c("InChIKey", "inchikey",
                                        "pubchem_inchikey"))
  smiles_col = .tanimoto_detect_col(raw, smiles_col,
                                    c("SMILES", "smiles", "CanonicalSMILES",
                                      "IsomericSMILES"))
  compound_id_col = .tanimoto_detect_col(raw, compound_id_col,
                                         c("compound_id", "CompoundID",
                                           "compound_key"))
  if (is.null(compound_col) && is.null(cid_col) && is.null(inchikey_col)) {
    stop("Could not detect a compound name, CID, or InChIKey column.",
         call. = FALSE)
  }

  group_cols = .tanimoto_validate_group_cols(raw, group_cols)
  compound_name = if (!is.null(compound_col)) {
    .uaf_squish_text(raw[[compound_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  cid = if (!is.null(cid_col)) .uaf_squish_text(raw[[cid_col]]) else {
    rep(NA_character_, nrow(raw))
  }
  inchikey = if (!is.null(inchikey_col)) {
    .uaf_squish_text(raw[[inchikey_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  smiles = if (!is.null(smiles_col)) {
    .uaf_squish_text(raw[[smiles_col]])
  } else {
    rep(NA_character_, nrow(raw))
  }
  compound_id = if (!is.null(compound_id_col)) {
    .uaf_squish_text(raw[[compound_id_col]])
  } else {
    .tanimoto_make_compound_id(inchikey, cid, compound_name)
  }
  missing_id = is.na(compound_id) | compound_id == ""
  compound_id[missing_id] = .tanimoto_make_compound_id(
    inchikey[missing_id], cid[missing_id], compound_name[missing_id]
  )

  normalized = data.frame(
    compound_id = compound_id,
    compound_name = compound_name,
    CID = cid,
    InChIKey = inchikey,
    SMILES = smiles,
    stringsAsFactors = FALSE
  )
  normalized = normalized[!is.na(normalized$compound_id) &
                            normalized$compound_id != "", , drop = FALSE]
  normalized = normalized[!duplicated(normalized$compound_id), , drop = FALSE]
  row.names(normalized) = NULL

  membership = .uaf_empty_table(c(group_cols, "group_id", "compound_id",
                                  "compound_name"))
  if (length(group_cols) > 0) {
    membership = raw[, group_cols, drop = FALSE]
    for (col in group_cols) membership[[col]] = .uaf_squish_text(membership[[col]])
    membership$group_id = .tanimoto_group_id(membership, group_cols)
    membership$compound_id = compound_id
    membership$compound_name = compound_name
    extra_cols = setdiff(names(raw), names(membership))
    for (col in extra_cols) membership[[col]] = raw[[col]]
    membership = membership[!is.na(membership$group_id) &
                              !is.na(membership$compound_id), , drop = FALSE]
    membership = membership[!duplicated(membership[, c("group_id",
                                                       "compound_id")]),
                            , drop = FALSE]
    row.names(membership) = NULL
  }
  list(compounds = normalized, membership = membership)
}

.tanimoto_detect_col = function(data, explicit, candidates) {
  if (!is.null(explicit)) {
    if (!explicit %in% names(data)) {
      stop("Column `", explicit, "` was not found.", call. = FALSE)
    }
    return(explicit)
  }
  hit = candidates[candidates %in% names(data)]
  if (length(hit) > 0) return(hit[[1]])
  lower = tolower(names(data))
  idx = match(tolower(candidates), lower, nomatch = 0)
  idx = idx[idx > 0]
  if (length(idx) > 0) names(data)[[idx[[1]]]] else NULL
}

.tanimoto_validate_group_cols = function(data, group_cols) {
  if (is.null(group_cols) || length(group_cols) < 1) return(character())
  missing = setdiff(group_cols, names(data))
  if (length(missing) > 0) {
    stop("Grouping columns were not found: ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  group_cols
}

.tanimoto_make_compound_id = function(inchikey, cid, name) {
  id = .tanimoto_sanitize_id(inchikey)
  missing = is.na(id) | id == ""
  cid_id = .tanimoto_sanitize_id(cid)
  id[missing & !is.na(cid_id)] = paste0("cid_", cid_id[missing & !is.na(cid_id)])
  missing = is.na(id) | id == ""
  name_id = .tanimoto_sanitize_id(name)
  id[missing & !is.na(name_id)] = paste0("name_", name_id[missing & !is.na(name_id)])
  id
}

.tanimoto_sanitize_id = function(x) {
  x = .uaf_squish_text(x)
  x = gsub("[^A-Za-z0-9]+", "_", x, perl = TRUE)
  x = gsub("^_+|_+$", "", x, perl = TRUE)
  x[x == ""] = NA_character_
  x
}

.tanimoto_group_id = function(data, group_cols) {
  if (length(group_cols) == 1) return(.uaf_squish_text(data[[group_cols]]))
  apply(data[, group_cols, drop = FALSE], 1, function(row) {
    paste(paste(group_cols, .uaf_squish_text(row), sep = "="),
          collapse = "|")
  })
}

.tanimoto_progress = function(progress_fun, phase, event, completed = NA_real_,
                              total = NA_real_, detail = NA_character_) {
  if (!is.function(progress_fun)) return(invisible(FALSE))
  payload = list(
    phase = phase,
    event = event,
    event_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    completed = completed,
    total = total,
    detail = detail
  )
  tryCatch(progress_fun(payload), error = function(error) {
    warning("Tanimoto progress callback failed: ", conditionMessage(error),
            call. = FALSE)
    invisible(FALSE)
  })
  invisible(TRUE)
}

.tanimoto_smoke_membership = function(membership, n = 25L,
                                      group_col = "species",
                                      compound_id_col = "compound_id") {
  membership = as.data.frame(membership, stringsAsFactors = FALSE)
  required = c(group_col, compound_id_col)
  missing = setdiff(required, names(membership))
  if (length(missing) > 0L) {
    stop("Smoke selection is missing columns: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  n = suppressWarnings(as.integer(n))
  if (length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be an integer of at least 2.", call. = FALSE)
  }
  groups = .uaf_squish_text(membership[[group_col]])
  compounds = .uaf_squish_text(membership[[compound_id_col]])
  keep = !is.na(groups) & groups != "" &
    !is.na(compounds) & compounds != ""
  membership = membership[keep, , drop = FALSE]
  groups = groups[keep]
  compounds = compounds[keep]
  if (length(unique(compounds)) < 2L) {
    stop("At least two unique structures are required for smoke selection.",
         call. = FALSE)
  }

  by_group = split(compounds, groups, drop = TRUE)
  by_group = lapply(by_group[order(names(by_group))], function(x) {
    sort(unique(x))
  })
  selected = character()
  max_rank = max(lengths(by_group))
  for (rank in seq_len(max_rank)) {
    for (group in names(by_group)) {
      candidates = by_group[[group]]
      if (rank > length(candidates)) next
      candidate = candidates[[rank]]
      if (!candidate %in% selected) selected = c(selected, candidate)
      if (length(selected) >= min(n, length(unique(compounds)))) break
    }
    if (length(selected) >= min(n, length(unique(compounds)))) break
  }
  selected = selected[seq_len(min(length(selected), n))]
  out = membership[
    .uaf_squish_text(membership[[compound_id_col]]) %in% selected,
    , drop = FALSE
  ]
  out$.smoke_selection_rank = match(
    .uaf_squish_text(out[[compound_id_col]]), selected
  )
  out = out[order(out$.smoke_selection_rank,
                  .uaf_squish_text(out[[group_col]])), , drop = FALSE]
  row.names(out) = NULL
  attr(out, "selected_compound_ids") = selected
  out
}

.tanimoto_resolve_cids = function(compounds, fetch, name_fallback,
                                  progress_fun = NULL,
                                  progress_every = 25L) {
  rows = lapply(seq_len(nrow(compounds)), function(i) {
    row = compounds[i, , drop = FALSE]
    cid = suppressWarnings(as.integer(row$CID))
    status = if (!is.na(cid)) "input_cid" else "not_queried"
    source = if (!is.na(cid)) "input_cid" else NA_character_
    url = NA_character_
    error = NA_character_
    if (is.na(cid) && !is.na(row$InChIKey) && row$InChIKey != "") {
      url = paste0(.pubchem_base_url(), "/pug/compound/inchikey/",
                   .pubchem_encode_path(row$InChIKey), "/cids/JSON")
      json = fetch(url)
      if (is.null(json)) {
        status = "inchikey_error"
        error = "PubChem InChIKey request failed."
      } else {
        cid = .tanimoto_json_cid(json)
        status = ifelse(is.na(cid), "inchikey_no_hit", "inchikey_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_inchikey")
      }
    }
    if (is.na(cid) && isTRUE(name_fallback) &&
        !is.na(row$compound_name) && row$compound_name != "") {
      url = paste0(.pubchem_base_url(), "/pug/compound/name/",
                   .pubchem_encode_path(row$compound_name), "/cids/JSON")
      json = fetch(url)
      if (is.null(json)) {
        status = "name_error"
        error = "PubChem name request failed."
      } else {
        cid = .tanimoto_json_cid(json)
        status = ifelse(is.na(cid), "name_no_hit", "name_resolved")
        source = ifelse(is.na(cid), NA_character_, "pubchem_name")
      }
    }
    out = data.frame(
      compound_id = row$compound_id,
      compound_name = row$compound_name,
      input_CID = row$CID,
      input_InChIKey = row$InChIKey,
      input_SMILES = row$SMILES,
      pubchem_cid = cid,
      cid_resolution_status = status,
      cid_resolution_source = source,
      cid_resolution_url = url,
      cid_resolution_error = error,
      stringsAsFactors = FALSE
    )
    if (i %% progress_every == 0L || i == nrow(compounds)) {
      .tanimoto_progress(
        progress_fun, "identity_resolution", "identities_completed",
        completed = i, total = nrow(compounds),
        detail = paste("Resolved", i, "of", nrow(compounds),
                       "submitted structures.")
      )
    }
    out
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.tanimoto_json_cid = function(json) {
  cid = tryCatch(json$IdentifierList$CID[[1]], error = function(error) NA)
  cid = suppressWarnings(as.integer(cid))
  if (length(cid) != 1 || is.na(cid)) NA_integer_ else cid
}

.tanimoto_fetch_fingerprints = function(resolution, fetch, chunk_size = 100,
                                        progress_fun = NULL,
                                        progress_every = 25L) {
  cids = unique(resolution$pubchem_cid[!is.na(resolution$pubchem_cid)])
  cols = c(
    "compound_id", "compound_name", "input_CID", "input_InChIKey",
    "input_SMILES", "pubchem_cid", "cid_resolution_status",
    "cid_resolution_source", "cid_resolution_url", "cid_resolution_error",
    "Title", "MolecularFormula", "InChIKey", "CanonicalSMILES",
    "IsomericSMILES", "Fingerprint2D", "fingerprint_status",
    "fingerprint_url"
  )
  if (length(cids) < 1) return(.uaf_empty_table(cols))
  chunks = split(cids, ceiling(seq_along(cids) / chunk_size))
  props = paste(c("Fingerprint2D", "CanonicalSMILES", "IsomericSMILES",
                  "InChIKey", "MolecularFormula", "Title"), collapse = ",")
  fetched = list()
  for (chunk_index in seq_along(chunks)) {
    chunk = chunks[[chunk_index]]
    url = paste0(.pubchem_base_url(), "/pug/compound/cid/",
                 paste(chunk, collapse = ","), "/property/", props, "/JSON")
    json = fetch(url)
    rows = tryCatch(json$PropertyTable$Properties,
                    error = function(error) list())
    if (length(rows) > 0) {
      fetched[[length(fetched) + 1L]] = do.call(rbind, lapply(rows, function(x) {
        data.frame(
          pubchem_cid = suppressWarnings(as.integer(.uaf_first_non_empty_text(x$CID))),
          Title = .uaf_first_non_empty_text(x$Title),
          MolecularFormula = .uaf_first_non_empty_text(x$MolecularFormula),
          InChIKey = .uaf_first_non_empty_text(x$InChIKey),
          CanonicalSMILES = .uaf_first_non_empty_text(x$CanonicalSMILES),
          IsomericSMILES = .uaf_first_non_empty_text(x$IsomericSMILES),
          Fingerprint2D = .uaf_first_non_empty_text(x$Fingerprint2D),
          fingerprint_status = ifelse(
            is.na(.uaf_first_non_empty_text(x$Fingerprint2D)),
            "fingerprint_missing", "fingerprint_resolved"
          ),
          fingerprint_url = url,
          stringsAsFactors = FALSE
        )
      }))
    }
    .tanimoto_progress(
      progress_fun, "fingerprint_resolution", "property_chunk_completed",
      completed = chunk_index, total = length(chunks),
      detail = paste("Completed PubChem fingerprint/property chunk",
                     chunk_index, "of", length(chunks), "for",
                     length(chunk), "CIDs.")
    )
  }
  if (length(fetched) < 1) return(.uaf_empty_table(cols))
  props = do.call(rbind, fetched)
  out = merge(resolution, props, by = "pubchem_cid", all.x = TRUE)
  out$fingerprint_status[is.na(out$fingerprint_status)] =
    "fingerprint_not_available"
  out = out[, intersect(cols, names(out)), drop = FALSE]
  row.names(out) = NULL
  out
}

.tanimoto_canonicalize_fingerprints = function(fingerprints, membership) {
  x = as.data.frame(fingerprints, stringsAsFactors = FALSE)
  required = c(
    "compound_id", "compound_name", "pubchem_cid", "input_CID",
    "input_InChIKey", "InChIKey", "Fingerprint2D",
    "cid_resolution_status", "cid_resolution_source", "fingerprint_status"
  )
  for (col in required) {
    if (!col %in% names(x)) x[[col]] = rep(NA_character_, nrow(x))
  }

  input_key = toupper(.uaf_squish_text(x$input_InChIKey))
  pubchem_key = toupper(.uaf_squish_text(x$InChIKey))
  has_input_key = !is.na(input_key) & input_key != ""
  has_pubchem_key = !is.na(pubchem_key) & pubchem_key != ""
  exact_key_match = has_input_key & has_pubchem_key &
    input_key == pubchem_key
  key_mismatch = has_input_key & has_pubchem_key &
    input_key != pubchem_key
  x$identity_match_status = ifelse(
    exact_key_match,
    "verified_input_inchikey",
    ifelse(
      key_mismatch,
      "input_pubchem_inchikey_mismatch",
      ifelse(
        has_input_key & !has_pubchem_key,
        "pubchem_inchikey_missing",
        ifelse(has_pubchem_key,
               "pubchem_inchikey_available_unverified_input",
               "no_inchikey_available")
      )
    )
  )
  x$canonical_compound_id = .tanimoto_make_compound_id(
    pubchem_key, x$pubchem_cid, x$compound_name
  )
  has_fingerprint = !is.na(x$Fingerprint2D) & x$Fingerprint2D != ""
  direct_cid_without_key = !has_input_key & !has_pubchem_key &
    tolower(.uaf_squish_text(x$cid_resolution_status)) == "input_cid"
  identity_consistent = exact_key_match |
    (!has_input_key & has_pubchem_key) | direct_cid_without_key
  initially_usable = has_fingerprint & identity_consistent &
    !is.na(x$canonical_compound_id) & x$canonical_compound_id != ""

  x$canonical_identity_conflict = FALSE
  candidate_keys = unique(x$canonical_compound_id[initially_usable])
  candidate_keys = candidate_keys[!is.na(candidate_keys) &
                                    candidate_keys != ""]
  for (key in candidate_keys) {
    idx = which(initially_usable & x$canonical_compound_id == key)
    cid_count = length(unique(.uaf_non_empty(x$pubchem_cid[idx])))
    inchikey_count = length(unique(.uaf_non_empty(pubchem_key[idx])))
    fingerprint_count = length(unique(.uaf_non_empty(x$Fingerprint2D[idx])))
    if (cid_count > 1L || inchikey_count > 1L || fingerprint_count > 1L) {
      x$canonical_identity_conflict[idx] = TRUE
    }
  }
  x$fingerprint_usable = initially_usable & !x$canonical_identity_conflict
  x$identity_exclusion_reason = ifelse(
    x$fingerprint_usable,
    NA_character_,
    ifelse(
      key_mismatch,
      "Input InChIKey does not match the InChIKey returned for the resolved PubChem CID.",
      ifelse(
        x$canonical_identity_conflict,
        "Canonical identity maps to conflicting PubChem CIDs, InChIKeys, or fingerprints.",
        ifelse(
          !has_fingerprint,
          "PubChem Fingerprint2D is unavailable.",
          ifelse(
            has_input_key & !has_pubchem_key,
            "PubChem did not return an InChIKey needed to verify the input identity.",
            ifelse(
              !identity_consistent,
              "Identity could not be verified from an input InChIKey or direct PubChem CID.",
              "Canonical compound identifier is unavailable."
            )
          )
        )
      )
    )
  )

  identity_map = data.frame(
    input_compound_id = x$compound_id,
    canonical_compound_id = x$canonical_compound_id,
    compound_name = x$compound_name,
    input_CID = x$input_CID,
    pubchem_cid = x$pubchem_cid,
    input_InChIKey = input_key,
    pubchem_InChIKey = pubchem_key,
    cid_resolution_status = x$cid_resolution_status,
    cid_resolution_source = x$cid_resolution_source,
    identity_match_status = x$identity_match_status,
    canonical_identity_conflict = x$canonical_identity_conflict,
    fingerprint_status = x$fingerprint_status,
    fingerprint_usable = x$fingerprint_usable,
    exclusion_reason = x$identity_exclusion_reason,
    stringsAsFactors = FALSE
  )

  usable = x[x$fingerprint_usable, , drop = FALSE]
  if (nrow(usable) > 0) {
    usable$input_compound_id = usable$compound_id
    usable = usable[order(usable$canonical_compound_id,
                          usable$input_compound_id,
                          usable$compound_name, na.last = TRUE), , drop = FALSE]
    groups = split(seq_len(nrow(usable)), usable$canonical_compound_id)
    source_ids = vapply(groups, function(idx) {
      .tanimoto_collapse_values(usable$input_compound_id[idx])
    }, character(1))
    aliases = vapply(groups, function(idx) {
      .tanimoto_collapse_values(usable$compound_name[idx])
    }, character(1))
    source_counts = vapply(groups, length, integer(1))
    compounds = usable[!duplicated(usable$canonical_compound_id), ,
                       drop = FALSE]
    group_index = match(compounds$canonical_compound_id, names(groups))
    compounds$source_compound_count = source_counts[group_index]
    compounds$source_compound_ids = source_ids[group_index]
    compounds$compound_aliases = aliases[group_index]
    compounds$compound_id = compounds$canonical_compound_id
    row.names(compounds) = NULL
  } else {
    compounds = x[FALSE, , drop = FALSE]
    compounds$source_compound_count = integer()
    compounds$source_compound_ids = character()
    compounds$compound_aliases = character()
  }

  membership_evidence = as.data.frame(membership, stringsAsFactors = FALSE)
  if (nrow(membership_evidence) > 0) {
    membership_evidence$tanimoto_input_compound_id =
      membership_evidence$compound_id
    map_index = match(membership_evidence$tanimoto_input_compound_id,
                      identity_map$input_compound_id)
    membership_evidence$canonical_compound_id =
      identity_map$canonical_compound_id[map_index]
    membership_evidence$identity_match_status =
      identity_map$identity_match_status[map_index]
    membership_evidence$fingerprint_usable =
      identity_map$fingerprint_usable[map_index]
    membership_evidence$identity_exclusion_reason =
      identity_map$exclusion_reason[map_index]
    membership_evidence$compound_id =
      membership_evidence$canonical_compound_id
    membership_evidence = membership_evidence[
      membership_evidence$fingerprint_usable %in% TRUE &
        !is.na(membership_evidence$compound_id), , drop = FALSE
    ]
    clean_membership = .tanimoto_collapse_membership(membership_evidence)
  } else {
    clean_membership = membership_evidence
  }

  excluded = identity_map[!(identity_map$fingerprint_usable %in% TRUE), ,
                          drop = FALSE]
  row.names(x) = row.names(identity_map) = row.names(excluded) = NULL
  list(
    fingerprints = x,
    compounds = compounds,
    identity_map = identity_map,
    excluded = excluded,
    membership_evidence = membership_evidence,
    membership = clean_membership
  )
}

.tanimoto_collapse_membership = function(membership) {
  if (!is.data.frame(membership) || nrow(membership) < 1) return(membership)
  if (!all(c("group_id", "compound_id") %in% names(membership))) {
    return(membership)
  }
  key = paste(membership$group_id, membership$compound_id, sep = "\r")
  groups = split(seq_len(nrow(membership)), key)
  rows = lapply(groups, function(idx) {
    part = membership[idx, , drop = FALSE]
    part = part[order(part$tanimoto_input_compound_id,
                      part$compound_name, na.last = TRUE), , drop = FALSE]
    row = part[1, , drop = FALSE]
    row$collapsed_membership_row_count = nrow(part)
    row$tanimoto_input_compound_ids = .tanimoto_collapse_values(
      part$tanimoto_input_compound_id
    )
    row$compound_aliases = .tanimoto_collapse_values(part$compound_name)
    row$source_record_count = if ("source_record_id" %in% names(part)) {
      length(.uaf_non_empty(unique(part$source_record_id)))
    } else {
      NA_integer_
    }
    row$source_record_ids = if ("source_record_id" %in% names(part)) {
      .tanimoto_collapse_values(part$source_record_id)
    } else {
      NA_character_
    }
    row$source_databases = if ("source_database" %in% names(part)) {
      .tanimoto_collapse_values(part$source_database)
    } else {
      NA_character_
    }
    row$evidence_tiers = if ("evidence_tier" %in% names(part)) {
      .tanimoto_collapse_values(part$evidence_tier)
    } else {
      NA_character_
    }
    row
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(out$group_id, out$compound_id), , drop = FALSE]
}

.tanimoto_collapse_values = function(x) {
  values = sort(unique(.uaf_non_empty(.uaf_squish_text(x))))
  if (length(values) < 1) return(NA_character_)
  paste(values, collapse = "; ")
}

.tanimoto_decode_fingerprints = function(fingerprints) {
  if (!requireNamespace("ChemmineR", quietly = TRUE)) {
    stop("ChemmineR is required to decode PubChem Fingerprint2D values.",
         call. = FALSE)
  }
  fp = fingerprints$Fingerprint2D
  names(fp) = fingerprints$compound_id
  keep = !is.na(fp) & fp != "" & !is.na(names(fp)) & names(fp) != ""
  fp = fp[keep]
  if (length(fp) < 2) {
    stop("Fewer than two fingerprints were available.", call. = FALSE)
  }
  bits = ChemmineR::fp2bit(fp, type = 2)
  storage.mode(bits) = "numeric"
  bits
}

.tanimoto_compound_pair_cols = function() {
  c("compound_pair_id", "compound_id_a", "compound_name_a",
    "pubchem_cid_a", "inchikey_a", "compound_id_b", "compound_name_b",
    "pubchem_cid_b", "inchikey_b", "tanimoto", "intersection_bits",
    "union_bits", "fingerprint_source")
}

.tanimoto_compound_pair_table = function(bits, compounds, block_size) {
  n = nrow(bits)
  counts = rowSums(bits)
  rows = list()
  starts = seq(1L, n - 1L, by = block_size)
  for (start in starts) {
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[block_idx, , drop = FALSE] %*% t(bits)
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      rows[[length(rows) + 1L]] =
        .tanimoto_compound_pair_rows(compounds, counts, inter[k, j],
                                     i = i, j = j)
    }
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.tanimoto_pair_shard_path = function(path, shard_index, multiple) {
  if (!isTRUE(multiple)) return(path)
  stem = sub("[.]csv[.]gz$", "", path, ignore.case = TRUE)
  if (identical(stem, path)) stem = sub("[.]gz$", "", path,
                                        ignore.case = TRUE)
  paste0(stem, "_part_", sprintf("%05d", shard_index), ".csv.gz")
}

.tanimoto_pair_writer = function(path, columns, expected_rows, shard_rows) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  expected_rows = suppressWarnings(as.numeric(expected_rows))
  if (length(expected_rows) != 1L || is.na(expected_rows) ||
      expected_rows < 0) {
    stop("`expected_rows` must be one non-negative number.", call. = FALSE)
  }
  shard_rows = suppressWarnings(as.numeric(shard_rows))
  if (length(shard_rows) != 1L || is.na(shard_rows) || shard_rows < 1) {
    stop("`shard_rows` must be a positive number or Inf.", call. = FALSE)
  }
  state = new.env(parent = emptyenv())
  state$multiple = is.finite(shard_rows) && expected_rows > shard_rows
  state$limit = if (state$multiple) floor(shard_rows) else Inf
  state$index = 0L
  state$con = NULL
  state$current_rows = 0
  state$partials = character()
  state$finals = character()
  state$rows = numeric()

  open_shard = function() {
    state$index = state$index + 1L
    final = .tanimoto_pair_shard_path(path, state$index, state$multiple)
    partial = paste0(final, ".partial")
    if (file.exists(partial)) unlink(partial, force = TRUE)
    state$con = gzfile(partial, "wt")
    writeLines(paste(columns, collapse = ","), state$con)
    state$current_rows = 0
    state$partials = c(state$partials, partial)
    state$finals = c(state$finals, final)
    state$rows = c(state$rows, 0)
    invisible(TRUE)
  }
  close_shard = function() {
    if (!is.null(state$con)) {
      close(state$con)
      state$con = NULL
    }
    invisible(TRUE)
  }
  write_rows = function(rows) {
    if (!is.data.frame(rows) || nrow(rows) < 1L) return(invisible(0L))
    offset = 1L
    while (offset <= nrow(rows)) {
      if (is.null(state$con)) open_shard()
      capacity = if (is.finite(state$limit)) {
        state$limit - state$current_rows
      } else {
        nrow(rows) - offset + 1L
      }
      take = min(nrow(rows) - offset + 1L, capacity)
      idx = seq.int(offset, length.out = take)
      utils::write.table(rows[idx, , drop = FALSE], state$con, sep = ",",
                         row.names = FALSE, col.names = FALSE, quote = TRUE,
                         na = "")
      state$current_rows = state$current_rows + take
      state$rows[[state$index]] = state$current_rows
      offset = offset + take
      if (is.finite(state$limit) && state$current_rows >= state$limit) {
        close_shard()
      }
    }
    invisible(nrow(rows))
  }
  abort = function() {
    close_shard()
    invisible(state$partials)
  }
  finalize = function() {
    if (state$index < 1L) open_shard()
    close_shard()
    if (sum(state$rows) != expected_rows) {
      stop("Pair writer produced ", sum(state$rows), " rows; expected ",
           expected_rows, ". Partial files were retained for audit.",
           call. = FALSE)
    }
    for (i in seq_along(state$partials)) {
      if (file.exists(state$finals[[i]])) unlink(state$finals[[i]],
                                                 force = TRUE)
      if (!file.rename(state$partials[[i]], state$finals[[i]])) {
        stop("Could not atomically publish pairwise output: ",
             state$finals[[i]], call. = FALSE)
      }
    }
    data.frame(
      Path = normalizePath(state$finals, winslash = "/", mustWork = FALSE),
      Rows = state$rows,
      Shard = seq_along(state$finals),
      ShardCount = length(state$finals),
      stringsAsFactors = FALSE
    )
  }
  list(write = write_rows, finalize = finalize, abort = abort)
}

.tanimoto_write_compound_pairs = function(bits, compounds, path, block_size,
                                           shard_rows = Inf,
                                           progress_fun = NULL,
                                           progress_every = 25L) {
  n = nrow(bits)
  counts = rowSums(bits)
  expected_rows = if (n >= 2L) choose(n, 2) else 0
  writer = .tanimoto_pair_writer(
    path, .tanimoto_compound_pair_cols(), expected_rows, shard_rows
  )
  completed = FALSE
  on.exit(if (!completed) writer$abort(), add = TRUE)
  if (n < 2L) {
    out = writer$finalize()
    completed = TRUE
    return(out)
  }
  starts = seq(1L, n - 1L, by = block_size)
  for (block_index in seq_along(starts)) {
    start = starts[[block_index]]
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[block_idx, , drop = FALSE] %*% t(bits)
    rows = lapply(seq_along(block_idx), function(k) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      .tanimoto_compound_pair_rows(compounds, counts, inter[k, j],
                                   i = i, j = j)
    })
    writer$write(do.call(rbind, rows))
    if (block_index %% progress_every == 0L ||
        block_index == length(starts)) {
      .tanimoto_progress(
        progress_fun, "compound_pairwise", "pair_blocks_completed",
        completed = block_index, total = length(starts),
        detail = paste("Completed compound-pair block", block_index, "of",
                       length(starts), ".")
      )
    }
  }
  out = writer$finalize()
  completed = TRUE
  out
}

.tanimoto_compound_pair_rows = function(compounds, counts, intersections,
                                        i, j) {
  intersection = as.numeric(intersections)
  union = counts[[i]] + counts[j] - intersection
  tanimoto = ifelse(union > 0, intersection / union, NA_real_)
  data.frame(
    compound_pair_id = paste(compounds$compound_id[[i]],
                             compounds$compound_id[j], sep = "__"),
    compound_id_a = compounds$compound_id[[i]],
    compound_name_a = compounds$compound_name[[i]],
    pubchem_cid_a = compounds$pubchem_cid[[i]],
    inchikey_a = compounds$InChIKey[[i]],
    compound_id_b = compounds$compound_id[j],
    compound_name_b = compounds$compound_name[j],
    pubchem_cid_b = compounds$pubchem_cid[j],
    inchikey_b = compounds$InChIKey[j],
    tanimoto = round(tanimoto, 6),
    intersection_bits = intersection,
    union_bits = union,
    fingerprint_source = "PubChem_Fingerprint2D",
    stringsAsFactors = FALSE
  )
}

.tanimoto_group_summary_cols = function(group_cols) {
  c("group_pair_id", "group_a", "group_b", "group_level",
    "compounds_a", "compounds_b", "compound_pair_count",
    "shared_compound_count", "mean_tanimoto", "median_tanimoto",
    "p95_tanimoto", "max_tanimoto", "top_compound_pair_ids",
    "top_compound_pairs", "fingerprint_source")
}

.tanimoto_group_pair_summary = function(bits, membership, compounds,
                                        group_cols, thresholds,
                                        top_n_pairs) {
  if (nrow(membership) < 2) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$fp_index) &
                            !is.na(membership$group_id), , drop = FALSE]
  groups = sort(unique(membership$group_id))
  if (length(groups) < 2) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  rows = list()
  for (i in seq_len(length(groups) - 1L)) {
    a = membership[membership$group_id == groups[[i]], , drop = FALSE]
    for (j in seq.int(i + 1L, length(groups))) {
      b = membership[membership$group_id == groups[[j]], , drop = FALSE]
      inter = bits[a$fp_index, , drop = FALSE] %*%
        t(bits[b$fp_index, , drop = FALSE])
      denom = outer(counts[a$fp_index], counts[b$fp_index], "+") - inter
      sim = ifelse(denom > 0, inter / denom, NA_real_)
      sim_vec = as.numeric(sim)
      sim_vec = sim_vec[!is.na(sim_vec)]
      if (length(sim_vec) < 1) next
      top = .tanimoto_top_pair_labels(sim, a, b, compound_names, top_n_pairs)
      thresholds_df = as.data.frame(stats::setNames(
        lapply(thresholds, function(th) sum(sim_vec >= th, na.rm = TRUE)),
        .tanimoto_threshold_cols(thresholds)
      ), stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] = cbind(data.frame(
        group_pair_id = paste(groups[[i]], groups[[j]], sep = "__"),
        group_a = groups[[i]],
        group_b = groups[[j]],
        group_level = paste(group_cols, collapse = "+"),
        compounds_a = nrow(a),
        compounds_b = nrow(b),
        compound_pair_count = length(sim_vec),
        shared_compound_count = length(intersect(a$compound_id,
                                                 b$compound_id)),
        mean_tanimoto = round(mean(sim_vec), 6),
        median_tanimoto = round(stats::median(sim_vec), 6),
        p95_tanimoto = round(unname(stats::quantile(sim_vec, 0.95,
                                                    na.rm = TRUE)), 6),
        max_tanimoto = round(max(sim_vec), 6),
        top_compound_pair_ids = paste(top$ids, collapse = "; "),
        top_compound_pairs = paste(top$labels, collapse = "; "),
        fingerprint_source = "PubChem_Fingerprint2D",
        stringsAsFactors = FALSE
      ), thresholds_df)
    }
  }
  if (length(rows) < 1) {
    return(.uaf_empty_table(c(.tanimoto_group_summary_cols(group_cols),
                              .tanimoto_threshold_cols(thresholds))))
  }
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out[order(-out$max_tanimoto, -out$mean_tanimoto, out$group_a,
            out$group_b), , drop = FALSE]
}

.tanimoto_comparable_pair_summaries = function(bits, membership, compounds,
                                                group_cols, thresholds,
                                                top_n_pairs, field,
                                                type = c("scope", "group")) {
  type = match.arg(type)
  cols = .comparable_pair_summary_cols(type, thresholds)
  if (!is.data.frame(membership) || nrow(membership) < 2L ||
      !field %in% names(membership) ||
      !"comparable_for_matrix" %in% names(membership)) {
    return(.uaf_empty_table(cols))
  }
  comparable = tolower(.uaf_squish_text(
    membership$comparable_for_matrix
  )) %in% c("yes", "true", "1")
  field_values = .uaf_squish_text(membership[[field]])
  comparable = comparable & !is.na(field_values) & field_values != "" &
    !tolower(field_values) %in% c("unknown", "mixed")
  values = sort(unique(field_values[comparable]))
  if (length(values) < 1L) return(.uaf_empty_table(cols))

  rows = lapply(values, function(value) {
    selected = membership[comparable & field_values == value, , drop = FALSE]
    summary = .tanimoto_group_pair_summary(
      bits, selected, compounds, group_cols, thresholds, top_n_pairs
    )
    if (nrow(summary) < 1L) return(NULL)
    support = lapply(summary$compound_pair_count, .comparable_support)
    out = data.frame(
      species_a = summary$group_a,
      species_b = summary$group_b,
      unordered_pair_key = paste(summary$group_a, summary$group_b,
                                 sep = " || "),
      comparison_type = type,
      comparison_value = value,
      compound_pair_count = summary$compound_pair_count,
      mean_tanimoto = summary$mean_tanimoto,
      median_tanimoto = summary$median_tanimoto,
      p95_tanimoto = summary$p95_tanimoto,
      max_tanimoto = summary$max_tanimoto,
      support_tier = vapply(support, `[[`, character(1), "tier"),
      low_support_caution = .bundle_yes_no(vapply(
        support, `[[`, logical(1), "low"
      )),
      support_note = vapply(support, `[[`, character(1), "note"),
      comparison_filter = paste0(field, " == '", value, "'"),
      caution = paste(
        "Filtered summary includes only compound pairs with matching",
        "comparable chemistry classification."
      ),
      stringsAsFactors = FALSE
    )
    threshold_cols = intersect(.tanimoto_threshold_cols(thresholds),
                               names(summary))
    for (col in threshold_cols) out[[col]] = summary[[col]]
    out[[field]] = value
    out
  })
  out = .plant_bind_tables(rows, cols)
  if (nrow(out) < 1L) return(.uaf_empty_table(cols))
  out = out[, cols, drop = FALSE]
  out[order(out$species_a, out$species_b, out$comparison_value), ,
      drop = FALSE]
}

.tanimoto_threshold_cols = function(thresholds) {
  paste0("compound_pair_count_ge_", gsub("\\.", "_", thresholds))
}

.tanimoto_top_pair_labels = function(sim, a, b, compound_names, top_n_pairs) {
  ord = order(as.numeric(sim), decreasing = TRUE, na.last = NA)
  if (length(ord) < 1) return(list(ids = character(), labels = character()))
  arr = arrayInd(utils::head(ord, top_n_pairs), .dim = dim(sim))
  ids = labels = character(nrow(arr))
  for (z in seq_len(nrow(arr))) {
    ca = a$compound_id[[arr[z, 1]]]
    cb = b$compound_id[[arr[z, 2]]]
    score = round(sim[arr[z, 1], arr[z, 2]], 3)
    ids[[z]] = paste0(ca, " | ", cb, "=", score)
    labels[[z]] = paste0(compound_names[[ca]], " | ",
                         compound_names[[cb]], "=", score)
  }
  list(ids = ids, labels = labels)
}

.tanimoto_group_pair_cols = function(group_cols) {
  c("group_compound_pair_id", "group_a", "compound_id_a",
    "compound_name_a", "group_b", "compound_id_b", "compound_name_b",
    "tanimoto", "intersection_bits", "union_bits", "fingerprint_source")
}

.tanimoto_group_pair_count = function(membership) {
  if (nrow(membership) < 2) return(0)
  total = choose(nrow(membership), 2)
  within = sum(vapply(split(membership$compound_id, membership$group_id),
                      function(x) choose(length(x), 2), numeric(1)))
  total - within
}

.tanimoto_group_compound_pair_table = function(bits, membership, compounds,
                                               block_size) {
  tf = tempfile(fileext = ".csv.gz")
  written = .tanimoto_write_group_compound_pairs(
    bits, membership, compounds, tf,
    block_size = block_size, shard_rows = Inf
  )
  tf = written$Path[[1]]
  utils::read.csv(gzfile(tf), stringsAsFactors = FALSE, check.names = FALSE)
}

.tanimoto_write_group_compound_pairs = function(bits, membership, compounds,
                                                path, block_size,
                                                shard_rows = Inf,
                                                progress_fun = NULL,
                                                progress_every = 25L) {
  membership$fp_index = match(membership$compound_id, rownames(bits))
  membership = membership[!is.na(membership$fp_index) &
                            !is.na(membership$group_id), , drop = FALSE]
  membership = membership[order(membership$group_id, membership$compound_id),
                          , drop = FALSE]
  counts = rowSums(bits)
  compound_names = stats::setNames(compounds$compound_name,
                                   compounds$compound_id)
  n = nrow(membership)
  expected_rows = .tanimoto_group_pair_count(membership)
  writer = .tanimoto_pair_writer(
    path, .tanimoto_group_pair_cols(character()), expected_rows, shard_rows
  )
  completed = FALSE
  on.exit(if (!completed) writer$abort(), add = TRUE)
  if (n < 2L) {
    out = writer$finalize()
    completed = TRUE
    return(out)
  }
  starts = seq(1L, n - 1L, by = block_size)
  for (block_index in seq_along(starts)) {
    start = starts[[block_index]]
    end = min(start + block_size - 1L, n - 1L)
    block_idx = start:end
    inter = bits[membership$fp_index[block_idx], , drop = FALSE] %*%
      t(bits[membership$fp_index, , drop = FALSE])
    rows = list()
    for (k in seq_along(block_idx)) {
      i = block_idx[[k]]
      j = seq.int(i + 1L, n)
      j = j[membership$group_id[j] != membership$group_id[[i]]]
      if (length(j) < 1) next
      intersection = as.numeric(inter[k, j])
      union = counts[membership$fp_index[[i]]] +
        counts[membership$fp_index[j]] - intersection
      tanimoto = ifelse(union > 0, intersection / union, NA_real_)
      ca = membership$compound_id[[i]]
      cb = membership$compound_id[j]
      rows[[length(rows) + 1L]] = data.frame(
        group_compound_pair_id = paste(
          paste(membership$group_id[[i]], ca, sep = "::"),
          paste(membership$group_id[j], cb, sep = "::"),
          sep = "__"
        ),
        group_a = membership$group_id[[i]],
        compound_id_a = ca,
        compound_name_a = compound_names[[ca]],
        group_b = membership$group_id[j],
        compound_id_b = cb,
        compound_name_b = unname(compound_names[cb]),
        tanimoto = round(tanimoto, 6),
        intersection_bits = intersection,
        union_bits = union,
        fingerprint_source = "PubChem_Fingerprint2D",
        stringsAsFactors = FALSE
      )
    }
    if (length(rows) > 0) {
      writer$write(do.call(rbind, rows))
    }
    if (block_index %% progress_every == 0L ||
        block_index == length(starts)) {
      .tanimoto_progress(
        progress_fun, "group_pairwise", "pair_blocks_completed",
        completed = block_index, total = length(starts),
        detail = paste("Completed cross-group pair block", block_index,
                       "of", length(starts), ".")
      )
    }
  }
  out = writer$finalize()
  completed = TRUE
  out
}

.tanimoto_plant_membership_table = function(plant_chemistry, level,
                                            include_review_required) {
  if (is.data.frame(plant_chemistry)) return(plant_chemistry)
  if (!is.list(plant_chemistry) ||
      !is.data.frame(plant_chemistry$PlantCompoundOccurrences) ||
      !is.data.frame(plant_chemistry$CompoundResolution)) {
    stop("`plant_chemistry` must be a plant phytochemistry result or a ",
         "plant-compound data frame.", call. = FALSE)
  }
  occ = plant_chemistry$PlantCompoundOccurrences
  cr = plant_chemistry$CompoundResolution
  for (col in c("compound_name_clean", "compound_name", "CID", "InChIKey",
                "SMILES", "MolecularFormula", "resolved",
                "resolution_source")) {
    if (!col %in% names(cr)) cr[[col]] = NA
  }
  cr$compound_id = .tanimoto_make_compound_id(cr$InChIKey, cr$CID,
                                              cr$compound_name_clean)
  cr = cr[cr$resolved %in% TRUE & !is.na(cr$compound_id), , drop = FALSE]
  if (isTRUE(include_review_required) &&
      is.data.frame(plant_chemistry$CompoundIdentityReview) &&
      nrow(plant_chemistry$CompoundIdentityReview) > 0) {
    review = plant_chemistry$CompoundIdentityReview
    if (!"compound_name_clean" %in% names(review)) {
      review$compound_name_clean = NA_character_
    }
    if (!"review_required" %in% names(review)) {
      review$review_required = TRUE
    }
    flags = stats::aggregate(
      list(identity_review_required = review$review_required %in% TRUE),
      by = list(compound_name_clean = .uaf_squish_text(review$compound_name_clean)),
      FUN = function(x) any(x %in% TRUE, na.rm = TRUE)
    )
    cr = merge(cr, flags, by = "compound_name_clean", all.x = TRUE)
    cr$identity_review_required[is.na(cr$identity_review_required)] = FALSE
  } else {
    cr$identity_review_required = FALSE
  }
  keep = c("compound_name_clean", "compound_id", "CID", "InChIKey", "SMILES",
           "MolecularFormula", "resolution_source",
           "identity_review_required")
  resolved = cr[, keep, drop = FALSE]
  names(resolved)[names(resolved) == "compound_id"] = "resolved_compound_id"
  joined = merge(occ, resolved, by = "compound_name_clean")
  if ("compound_id" %in% names(joined)) {
    joined$source_compound_id = joined$compound_id
  }
  joined$compound_id = joined$resolved_compound_id
  if (!"compound_name" %in% names(joined) && "compound_name.x" %in% names(joined)) {
    joined$compound_name = joined$compound_name.x
  }
  if (!level %in% names(joined)) {
    stop("Plant occurrence table does not contain `", level, "`.",
         call. = FALSE)
  }
  joined = joined[!is.na(joined[[level]]) & !is.na(joined$compound_id), ,
                  drop = FALSE]
  joined = joined[!duplicated(joined[, c(level, "compound_id")]), ,
                  drop = FALSE]
  row.names(joined) = NULL
  joined
}

.tanimoto_label_plant_summary = function(summary, level) {
  if (!is.data.frame(summary) || nrow(summary) < 1) return(summary)
  names(summary)[names(summary) == "group_a"] = paste0(level, "_a")
  names(summary)[names(summary) == "group_b"] = paste0(level, "_b")
  names(summary)[names(summary) == "group_pair_id"] =
    paste0(level, "_pair_id")
  summary
}

.tanimoto_label_plant_pairs = function(pairs, level) {
  if (!is.data.frame(pairs) || nrow(pairs) < 1) return(pairs)
  names(pairs)[names(pairs) == "group_a"] = paste0(level, "_a")
  names(pairs)[names(pairs) == "group_b"] = paste0(level, "_b")
  names(pairs)[names(pairs) == "group_compound_pair_id"] =
    paste0(level, "_compound_pair_id")
  pairs
}

.tanimoto_empty_manifest = function() {
  .uaf_empty_table(c("Table", "Path", "Rows", "Description"))
}

.tanimoto_manifest_row = function(table, path, rows, description) {
  data.frame(Table = table,
             Path = normalizePath(path, winslash = "/", mustWork = FALSE),
             Rows = rows,
             Description = description,
             stringsAsFactors = FALSE)
}

.tanimoto_manifest_rows = function(table, written, description) {
  if (!is.data.frame(written) ||
      !all(c("Path", "Rows") %in% names(written))) {
    stop("Pair writer did not return a valid output manifest.",
         call. = FALSE)
  }
  data.frame(
    Table = rep(table, nrow(written)),
    Path = normalizePath(written$Path, winslash = "/", mustWork = FALSE),
    Rows = written$Rows,
    Description = rep(description, nrow(written)),
    stringsAsFactors = FALSE
  )
}
