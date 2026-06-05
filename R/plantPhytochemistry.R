#' Species-first plant phytochemistry resolver
#'
#' @description
#' `resolvePlantPhytochemistry()` builds a species-first phytochemistry result
#' from plant names, public-provider adapters, literature candidates, and
#' optional curated species-compound rows. The resolver is intentionally
#' conservative: no compound is fabricated, literature co-mentions remain
#' candidate evidence, and provider no-hit states are reported through
#' diagnostics.
#'
#' @param plants Character vector or data frame of plant names. If a data frame
#' is supplied, a `species` column is preferred; otherwise the first column is
#' treated as the plant name. Optional `genus` and `family` columns are used
#' when present.
#' @param sources Character vector of provider names. Supported names are
#' `"lotus"`, `"knapsack"`, `"npass"`, `"pubchem"`, `"pubmed"`, and
#' `"pubtator"`.
#' @param taxon_fallback Character vector of fallback ranks to record in query
#' planning. Species-level evidence is queried first. Providers that support
#' broader lookup may also use genus or family terms; fallback evidence is
#' labeled separately from direct species evidence.
#' @param enrich_compounds Logical. If `TRUE`, attempt compound enrichment after
#' occurrence discovery.
#' @param chemical_library Optional chemical library passed to `categorate()`
#' when compound enrichment is requested. If omitted, compound enrichment falls
#' back to `pubchemProfile()` and returns PubChem-derived categorate-like
#' analysis tables without FMCS library matching.
#' @param detail One of `"research"`, `"full"`, or `"none"`.
#' @param cache Logical. If `TRUE`, provider requests may use cache files.
#' @param cache_dir Cache directory. Defaults to a uafR cache directory when
#' needed.
#' @param throttle Seconds to wait between uncached provider requests.
#' @param ncbi_email Optional NCBI email. Defaults to `Sys.getenv("NCBI_EMAIL")`.
#' @param ncbi_tool Optional NCBI tool name. Defaults to `Sys.getenv("NCBI_TOOL",
#' "uafR")`.
#' @param ncbi_api_key Optional NCBI API key. Defaults to
#' `Sys.getenv("NCBI_API_KEY")`. It is used in request URLs but is not persisted
#' in result tables.
#' @param max_pubmed_records Maximum PubMed records per plant to request.
#' @param max_provider_records Maximum occurrence rows to retain per plant from
#' non-PubMed provider adapters.
#' @param min_confidence Minimum confidence for summary/matrix features.
#' @param enrichment_batch_size Maximum number of unique compounds per
#' PubChem-only enrichment batch. Use `Inf` for one batch.
#' @param resume_enrichment Logical. If `TRUE`, reuse saved PubChem-only
#' enrichment batch files when available.
#' @param progress Logical. If `TRUE`, print simple enrichment progress
#' messages.
#' @param strict Logical. If `TRUE`, validation treats optional missing fields
#' more strictly.
#' @param curated_data Optional local species-compound table to standardize and
#' include.
#' @param provider_results Optional named list of mocked or pre-fetched provider
#' results. This is intended for tests, offline workflows, and future cached
#' provider download parsers.
#' @param enrichment_fun Optional function used instead of `categorate()` for
#' compound enrichment. It should accept a character vector of compounds and
#' return a categorate-like list.
#' @param request_fun Optional URL request function used by live-capable
#' provider adapters in tests or advanced use.
#' @param pubtator_request_fun Optional URL request function used by the
#' PubTator adapter in tests or advanced use.
#' @param refresh Logical. Reserved for future provider cache refresh support.
#' @param ... Additional arguments passed to `categorate()`, `enrichment_fun`,
#' or the PubChem-only enrichment fallback. Use `compound_request_fun` in `...`
#' to mock or customize compound-centered PubChem requests.
#'
#' @return A list with class `"uaf_plant_phytochemistry"` containing normalized
#' query, occurrence, literature, enrichment, summary, matrix, validation,
#' dictionary, and provenance tables.
#'
#' @examples
#' \dontrun{
#' plants = c("Salix nigra", "Camellia sinensis", "Zea mays")
#' phyto = resolvePlantPhytochemistry(
#'   plants = plants,
#'   sources = c("lotus", "pubmed", "pubtator"),
#'   taxon_fallback = c("species", "genus"),
#'   enrich_compounds = FALSE
#' )
#' phyto$SpeciesChemistrySummary
#' validatePlantPhytochemistryResult(phyto)$Summary
#' }
#'
#' @export
resolvePlantPhytochemistry = function(plants,
                                      sources = c("lotus", "knapsack",
                                                  "npass", "pubchem",
                                                  "pubmed", "pubtator"),
                                      taxon_fallback = c("species", "genus"),
                                      enrich_compounds = TRUE,
                                      chemical_library = NULL,
                                      detail = c("research", "full", "none"),
                                      cache = TRUE,
                                      cache_dir = NULL,
                                      throttle = 0.2,
                                      ncbi_email = Sys.getenv("NCBI_EMAIL", ""),
                                      ncbi_tool = Sys.getenv("NCBI_TOOL", "uafR"),
                                      ncbi_api_key = Sys.getenv("NCBI_API_KEY", ""),
                                      max_pubmed_records = 50,
                                      max_provider_records = max_pubmed_records,
                                      min_confidence = "medium",
                                      enrichment_batch_size = Inf,
                                      resume_enrichment = TRUE,
                                      progress = interactive(),
                                      strict = FALSE,
                                      curated_data = NULL,
                                      provider_results = NULL,
                                      enrichment_fun = NULL,
                                      request_fun = NULL,
                                      pubtator_request_fun = NULL,
                                      refresh = FALSE,
                                      ...) {
  detail = match.arg(detail)
  plant_queries = .plant_queries(plants, taxon_fallback)
  plant_resolution = .plant_name_resolution(plant_queries)
  sources = .plant_normalize_sources(sources)
  cache_dir = .plant_cache_dir(cache_dir)

  dispatch = .plant_provider_dispatch(
    plant_queries = plant_queries,
    sources = sources,
    cache = cache,
    cache_dir = cache_dir,
    throttle = throttle,
    ncbi_email = ncbi_email,
    ncbi_tool = ncbi_tool,
    ncbi_api_key = ncbi_api_key,
    max_pubmed_records = max_pubmed_records,
    max_provider_records = max_provider_records,
    provider_results = provider_results,
    request_fun = request_fun,
    pubtator_request_fun = pubtator_request_fun,
    refresh = refresh
  )

  occurrence_parts = list(dispatch$PlantCompoundOccurrences)
  if (!is.null(curated_data)) {
    occurrence_parts[[length(occurrence_parts) + 1]] =
      standardizePlantCompoundIntake(curated_data)
  }
  occurrences = .plant_bind_occurrences(occurrence_parts)
  occurrences = .plant_match_occurrences_to_queries(occurrences, plant_queries)

  enrichment = if (isTRUE(enrich_compounds) && detail != "none") {
    enrichPlantCompounds(occurrences,
                         chemical_library = chemical_library,
                         detail = detail,
                         cache = cache,
                         cache_dir = file.path(cache_dir, "compound_enrichment"),
                         throttle = throttle,
                         enrichment_fun = enrichment_fun,
                         batch_size = enrichment_batch_size,
                         resume = resume_enrichment,
                         progress = progress,
                         ...)
  } else {
    list(CategorateResult = NULL,
         CompoundResolution = .plant_compound_resolution(occurrences, NULL),
         TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
         Provenance = .plant_provenance("compound_enrichment", "not_run",
                                        NA_character_, NA_character_, 0,
                                        "Compound enrichment was not requested."))
  }

  summary = summarizePlantPhytochemistry(
    plant_compounds = occurrences,
    categorate_result = enrichment$CategorateResult,
    compound_resolution = enrichment$CompoundResolution,
    plant_queries = plant_queries,
    provider_diagnostics = dispatch$ProviderDiagnostics
  )
  matrix = plantPhytochemistryMatrix(
    list(PlantCompoundOccurrences = occurrences,
         CategorateResult = enrichment$CategorateResult),
    level = "species",
    profile = "core",
    mode = "binary",
    min_confidence = min_confidence
  )

  provenance = .plant_bind_tables(list(
    dispatch$Provenance,
    enrichment$Provenance,
    .plant_provenance("resolvePlantPhytochemistry", "uafR",
                      paste(plant_queries$query_plant, collapse = "; "),
                      NA_character_, nrow(occurrences),
                      "Species-first plant phytochemistry result assembled.")
  ), .plant_provenance_cols())

  out = list(
    PlantQueries = plant_queries,
    PlantNameResolution = plant_resolution,
    ProviderDiagnostics = dispatch$ProviderDiagnostics,
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = dispatch$LiteratureCandidates,
    CompoundResolution = enrichment$CompoundResolution,
    CategorateResult = enrichment$CategorateResult,
    SpeciesChemistrySummary = summary,
    SpeciesChemistryMatrix = matrix,
    TraitEvidence = enrichment$TraitEvidence,
    Validation = list(Summary = data.frame(),
                      TableQuality = data.frame(),
                      ProviderDiagnostics = data.frame(),
                      Issues = data.frame(),
                      DataDictionary = data.frame()),
    DataDictionary = plantPhytochemistrySchema(),
    Provenance = provenance
  )
  class(out) = c("uaf_plant_phytochemistry", class(out))
  validation = validatePlantPhytochemistryResult(out, strict = strict)
  out$Validation = validation
  out
}

#' Describe plant phytochemistry output schemas
#'
#' @param tables Optional table names to return. If `NULL`, all schemas are
#' returned.
#'
#' @return A data frame describing normalized table, column, type, required,
#' role, allowed-value, and description contracts.
#'
#' @export
plantPhytochemistrySchema = function(tables = NULL) {
  dictionary = .plant_data_dictionary()
  tables = .uaf_non_empty(tables)
  if (length(tables) > 0) {
    dictionary = dictionary[dictionary$Table %in% tables, , drop = FALSE]
  }
  row.names(dictionary) = NULL
  dictionary
}

#' Validate plant phytochemistry results
#'
#' @param x A result returned by `resolvePlantPhytochemistry()` or a named list
#' containing plant phytochemistry tables.
#' @param strict Logical. If `TRUE`, missing optional documented columns are
#' reported as warnings.
#'
#' @return A list with `Summary`, `TableQuality`, `ProviderDiagnostics`,
#' `Issues`, and `DataDictionary`.
#'
#' @export
validatePlantPhytochemistryResult = function(x, strict = FALSE) {
  dictionary = plantPhytochemistrySchema()
  expected_tables = unique(dictionary$Table)
  issue_rows = list()
  quality_rows = list()

  if (!is.list(x) || is.data.frame(x)) {
    issues = .plant_validation_issue("error", NA_character_, NA_character_,
                                     "Input is not a plant phytochemistry result list",
                                     "list", class(x)[[1]], NA_integer_,
                                     NA_character_)
    return(.plant_validation_result(.uaf_empty_table(.plant_table_quality_cols()),
                                    .uaf_empty_table(.plant_provider_diagnostic_cols()),
                                    issues, dictionary))
  }

  for (table_name in expected_tables) {
    table_dictionary = dictionary[dictionary$Table == table_name, ,
                                  drop = FALSE]
    table = x[[table_name]]
    quality_rows[[length(quality_rows) + 1]] =
      .plant_table_quality(table_name, table, table_dictionary, strict)
    issues = .plant_validate_table(table_name, table, table_dictionary, strict)
    if (nrow(issues) > 0) issue_rows[[length(issue_rows) + 1]] = issues
  }

  issues = .plant_bind_tables(issue_rows, .plant_validation_issue_cols())
  issues = .plant_extra_validation_issues(x, issues)
  table_quality = .plant_bind_tables(quality_rows, .plant_table_quality_cols())
  provider_diagnostics = if (is.data.frame(x$ProviderDiagnostics)) {
    x$ProviderDiagnostics
  } else {
    .uaf_empty_table(.plant_provider_diagnostic_cols())
  }

  .plant_validation_result(table_quality, provider_diagnostics, issues,
                           dictionary)
}

#' Standardize curated plant-compound intake
#'
#' @param x Data frame containing at least `species` and `compound_name`.
#' Recommended columns include `source_database`, `citation_or_url`, and
#' `evidence_tier`.
#'
#' @return A normalized `PlantCompoundOccurrences` data frame.
#'
#' @export
standardizePlantCompoundIntake = function(x) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data frame with plant-compound rows.", call. = FALSE)
  }
  names(x) = .plant_normalize_column_names(names(x))
  if (!"species" %in% names(x) || !"compound_name" %in% names(x)) {
    stop("Curated plant-compound intake requires `species` and ",
         "`compound_name` columns.", call. = FALSE)
  }

  out = .plant_empty_occurrences()
  if (nrow(x) < 1) return(out)
  species = .uaf_squish_text(x$species)
  genus = if ("genus" %in% names(x)) .uaf_squish_text(x$genus) else
    .plant_genus(species)
  family = if ("family" %in% names(x)) .uaf_squish_text(x$family) else
    rep(NA_character_, nrow(x))
  source_database = .plant_col_or_default(x, "source_database", "manual")
  citation = .plant_col_or_default(x, "citation_or_url", NA_character_)
  evidence_tier = .plant_col_or_default(x, "evidence_tier", "manual_curated")
  confidence = .plant_evidence_tier_confidence(evidence_tier)

  out = data.frame(
    query_plant = species,
    query_plant_clean = .plant_clean_name(species),
    matched_taxon = species,
    matched_rank = "species",
    species = species,
    genus = genus,
    family = family,
    compound_name = .uaf_squish_text(x$compound_name),
    compound_name_clean = .plant_clean_compound(x$compound_name),
    compound_id = .plant_col_or_default(x, "compound_id", NA_character_),
    compound_id_type = .plant_col_or_default(x, "compound_id_type", NA_character_),
    source_database = source_database,
    source_record_id = .plant_col_or_default(x, "source_record_id", NA_character_),
    evidence_text = .plant_col_or_default(x, "evidence_note", NA_character_),
    evidence_url = citation,
    reference_id = .plant_col_or_default(x, "reference_id", NA_character_),
    pmid = .plant_col_or_default(x, "pmid", NA_character_),
    doi = .plant_col_or_default(x, "doi", NA_character_),
    plant_part = .plant_col_or_default(x, "plant_part", NA_character_),
    tissue = .plant_col_or_default(x, "tissue", NA_character_),
    method = .plant_col_or_default(x, "method", NA_character_),
    occurrence_type = .plant_col_or_default(x, "occurrence_type", "reported"),
    retrieved_at = .plant_timestamp(),
    confidence = confidence,
    curation_flag = .plant_col_or_default(x, "curation_status", "curated"),
    evidence_tier = evidence_tier,
    stringsAsFactors = FALSE
  )
  out = .plant_normalize_occurrences(out, source_hint = "manual")
  row.names(out) = NULL
  out
}

#' Enrich plant compounds with existing uafR compound workflows
#'
#' @param plant_compounds A normalized or curated plant-compound table.
#' @param chemical_library Optional library passed to `categorate()`.
#' @param detail One of `"research"`, `"full"`, or `"none"`.
#' @param cache Logical.
#' @param cache_dir Cache directory.
#' @param throttle Request throttle.
#' @param enrichment_fun Optional injected enrichment function for tests or
#' cached workflows.
#' @param pubchem_fun Optional replacement for `pubchemProfile()` used by the
#' PubChem-only fallback. Intended for tests and advanced cached workflows.
#' @param compound_request_fun Optional request function passed to
#' `pubchemProfile()` when the PubChem-only fallback is used.
#' @param batch_size Maximum number of unique compounds per PubChem-only
#' enrichment batch. Use `Inf` for one batch.
#' @param resume Logical. If `TRUE`, reuse saved PubChem-only batch `.rds`
#' files in `cache_dir`.
#' @param progress Logical. If `TRUE`, print simple progress messages during
#' PubChem-only enrichment.
#' @param ... Additional arguments passed to `categorate()`, `enrichment_fun`,
#' or the PubChem-only enrichment fallback.
#'
#' @return A list with `CategorateResult`, `CompoundResolution`,
#' `TraitEvidence`, and `Provenance`.
#'
#' @export
enrichPlantCompounds = function(plant_compounds,
                                chemical_library = NULL,
                                detail = c("research", "full", "none"),
                                cache = TRUE,
                                cache_dir = NULL,
                                throttle = 0.2,
                                enrichment_fun = NULL,
                                pubchem_fun = NULL,
                                compound_request_fun = NULL,
                                batch_size = Inf,
                                resume = TRUE,
                                progress = interactive(),
                                ...) {
  detail = match.arg(detail)
  occurrences = if (is.data.frame(plant_compounds) &&
                    all(.plant_occurrence_cols() %in% names(plant_compounds))) {
    .plant_normalize_occurrences(plant_compounds)
  } else {
    standardizePlantCompoundIntake(plant_compounds)
  }
  compounds = unique(.uaf_non_empty(occurrences$compound_name))
  if (length(compounds) < 1 || detail == "none") {
    return(list(
      CategorateResult = NULL,
      CompoundResolution = .plant_compound_resolution(occurrences, NULL),
      TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
      Provenance = .plant_provenance("compound_enrichment", "not_run",
                                     paste(compounds, collapse = "; "),
                                     NA_character_, 0,
                                     "No compounds were available for enrichment or detail = 'none'.")
    ))
  }

  categorate_result = tryCatch({
    if (!is.null(enrichment_fun)) {
      enrichment_fun(compounds, ...)
    } else if (!is.null(chemical_library)) {
      categorate(compounds = compounds,
                 chemical_library = chemical_library,
                 detail = detail,
                 cache = cache,
                 cache_dir = cache_dir,
                 throttle = throttle,
                 ...)
    } else {
      .plant_pubchem_only_enrichment(
        compounds = compounds,
        detail = detail,
        cache = cache,
        cache_dir = cache_dir,
        throttle = throttle,
        pubchem_fun = pubchem_fun,
        request_fun = compound_request_fun,
        batch_size = batch_size,
        resume = resume,
        progress = progress,
        ...
      )
    }
  }, error = function(error) {
    warning("Plant compound enrichment failed: ", conditionMessage(error),
            call. = FALSE)
    NULL
  })

  list(
    CategorateResult = categorate_result,
    CompoundResolution = .plant_compound_resolution(occurrences,
                                                    categorate_result),
    TraitEvidence = .plant_trait_evidence(occurrences, categorate_result),
    Provenance = .plant_provenance("compound_enrichment",
                                   .plant_enrichment_source(categorate_result),
                                   paste(compounds, collapse = "; "),
                                   NA_character_, length(compounds),
                                   ifelse(is.null(categorate_result),
                                          "Compound enrichment not available.",
                                          .plant_enrichment_note(categorate_result)))
  )
}

#' Summarize plant phytochemistry
#'
#' @param plant_compounds Plant-compound occurrence table or a
#' `uaf_plant_phytochemistry` result.
#' @param categorate_result Optional categorate-like enrichment result.
#' @param compound_resolution Optional compound resolution table.
#' @param plant_queries Optional plant query table.
#' @param provider_diagnostics Optional provider diagnostics table.
#'
#' @return Species-level summary table.
#'
#' @export
summarizePlantPhytochemistry = function(plant_compounds,
                                        categorate_result = NULL,
                                        compound_resolution = NULL,
                                        plant_queries = NULL,
                                        provider_diagnostics = NULL) {
  if (inherits(plant_compounds, "uaf_plant_phytochemistry")) {
    x = plant_compounds
    plant_compounds = x$PlantCompoundOccurrences
    categorate_result = x$CategorateResult
    compound_resolution = x$CompoundResolution
    plant_queries = x$PlantQueries
    provider_diagnostics = x$ProviderDiagnostics
  }
  occurrences = .plant_normalize_occurrences(plant_compounds)
  if (is.null(plant_queries)) plant_queries = .plant_queries(unique(occurrences$species))
  if (is.null(compound_resolution)) {
    compound_resolution = .plant_compound_resolution(occurrences,
                                                     categorate_result)
  }

  base_species = unique(.uaf_non_empty(c(plant_queries$species,
                                         occurrences$species)))
  if (length(base_species) < 1) return(.uaf_empty_table(.plant_summary_cols()))
  rows = lapply(base_species, function(species) {
    group = occurrences[occurrences$species == species, , drop = FALSE]
    query = plant_queries[plant_queries$species == species, , drop = FALSE]
    resolved = compound_resolution[
      compound_resolution$compound_name_clean %in% group$compound_name_clean,
      , drop = FALSE
    ]
    source_count = length(unique(.uaf_non_empty(group$source_database)))
    direct_species = group[group$matched_rank == "species", , drop = FALSE]
    genus_fallback = group[group$matched_rank == "genus", , drop = FALSE]
    family_fallback = group[group$matched_rank == "family", , drop = FALSE]
    literature = group[grepl("literature|pubtator", group$evidence_tier,
                             ignore.case = TRUE), , drop = FALSE]
    database_rows = group[!grepl("literature|pubtator", group$evidence_tier,
                                 ignore.case = TRUE), , drop = FALSE]
    ready = group[group$analysis_ready == "Yes", , drop = FALSE]
    quality = suppressWarnings(as.numeric(group$evidence_quality_score))
    context_known = group[group$biological_context_status !=
                            "context_missing", , drop = FALSE]
    trait_info = .plant_summary_trait_info(group, categorate_result)
    data.frame(
      species = species,
      species_slug = .plant_slug(species),
      genus = .uaf_first_non_empty_text(query$genus, .plant_genus(species)),
      family = .uaf_first_non_empty_text(query$family, group$family),
      query_status = ifelse(nrow(group) > 0, "records_found", "no_public_records"),
      compound_count = length(unique(.uaf_non_empty(group$compound_name_clean))),
      resolved_compound_count = sum(resolved$resolved %in% TRUE),
      unresolved_compound_count = sum(!(resolved$resolved %in% TRUE)),
      analysis_ready_compound_count = length(unique(.uaf_non_empty(ready$compound_name_clean))),
      direct_species_compound_count = length(unique(.uaf_non_empty(direct_species$compound_name_clean))),
      genus_level_compound_count = length(unique(.uaf_non_empty(genus_fallback$compound_name_clean))),
      family_level_compound_count = length(unique(.uaf_non_empty(family_fallback$compound_name_clean))),
      literature_candidate_count = nrow(literature),
      database_occurrence_count = nrow(database_rows),
      direct_reported_occurrence_count = sum(group$occurrence_status ==
                                               "direct_reported"),
      curated_reported_occurrence_count = sum(group$occurrence_status ==
                                                "curated_reported"),
      candidate_occurrence_count = sum(group$occurrence_status == "candidate"),
      taxon_fallback_occurrence_count = sum(group$occurrence_status ==
                                              "taxon_fallback"),
      context_known_occurrence_count = nrow(context_known),
      mean_evidence_quality_score = ifelse(length(quality[is.finite(quality)]) > 0,
                                           round(mean(quality[is.finite(quality)]), 3),
                                           NA_real_),
      analysis_ready_fraction = ifelse(nrow(group) > 0,
                                       round(nrow(ready) / nrow(group), 3), 0),
      plant_part_groups = .pubchem_collapse(group$plant_part_group[
        !group$plant_part_group %in% c("unknown", "extract_unspecified")
      ]),
      tissue_groups = .pubchem_collapse(group$tissue_group[
        !group$tissue_group %in% c("unknown", "extract_unspecified")
      ]),
      method_groups = .pubchem_collapse(group$method_group[
        !group$method_group %in% c("unknown")
      ]),
      evidence_tier_summary = .plant_count_summary(group$evidence_tier),
      source_database_count = source_count,
      source_databases = .pubchem_collapse(group$source_database),
      natural_product_superclasses = trait_info$superclasses,
      natural_product_classes = trait_info$classes,
      natural_product_subclasses = trait_info$subclasses,
      dominant_compound_classes = trait_info$dominant_classes,
      kingdoms_observed = trait_info$kingdoms,
      families_observed = .pubchem_collapse(unique(.uaf_non_empty(group$family))),
      is_plant_occurring = any(group$occurrence_status %in%
                                 c("direct_reported", "curated_reported")),
      volatile_proxy_fraction = trait_info$volatile_fraction,
      lipophilic_fraction = trait_info$lipophilic_fraction,
      oxygenated_fraction = trait_info$oxygenated_fraction,
      nitrogenous_fraction = trait_info$nitrogenous_fraction,
      sulfur_containing_fraction = trait_info$sulfur_fraction,
      halogenated_fraction = trait_info$halogen_fraction,
      kegg_pathway_groups = trait_info$kegg_pathway_groups,
      chemical_trait_count = trait_info$chemical_trait_count,
      high_confidence_trait_count = trait_info$high_confidence_trait_count,
      source_coverage_score = .plant_source_coverage_score(source_count,
                                                           provider_diagnostics),
      uafR_validation_status = .plant_categorate_validation_status(categorate_result),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

#' Build plant phytochemistry matrices
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param level One of `"species"`, `"genus"`, or `"family"`.
#' @param profile One of `"core"`, `"full"`, `"ecology"`, `"metabolism"`,
#' `"safety"`, or `"bioactivity"`.
#' @param mode One of `"binary"`, `"count"`, or `"confidence"`.
#' @param min_confidence Minimum confidence to include.
#' @param max_traits Maximum number of columns to retain.
#'
#' @return Wide matrix-like data frame.
#'
#' @export
plantPhytochemistryMatrix = function(x,
                                     level = c("species", "genus", "family"),
                                     profile = c("core", "full", "ecology",
                                                 "metabolism", "safety",
                                                 "bioactivity"),
                                     mode = c("binary", "count", "confidence"),
                                     min_confidence = "medium",
                                     max_traits = Inf) {
  level = match.arg(level)
  profile = match.arg(profile)
  mode = match.arg(mode)
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else if (is.list(x) && is.data.frame(x$PlantCompoundOccurrences)) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  categorate_result = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$CategorateResult
  } else if (is.list(x)) {
    x$CategorateResult
  } else {
    NULL
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(c(level)))
  }
  keep = .plant_confidence_score(occurrences$confidence) >=
    .plant_confidence_score(min_confidence)
  occurrences = occurrences[keep, , drop = FALSE]
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(c(level)))
  }
  id = occurrences[[level]]
  id[is.na(id) | id == ""] = occurrences$species[is.na(id) | id == ""]
  long = data.frame(
    id = id,
    trait = paste0("compound__", .plant_matrix_key(occurrences$compound_name_clean)),
    value = .plant_confidence_score(occurrences$confidence),
    stringsAsFactors = FALSE
  )
  if (profile %in% c("core", "full", "ecology")) {
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("source__", .plant_matrix_key(occurrences$source_database)),
      value = .plant_confidence_score(occurrences$confidence),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("evidence__", .plant_matrix_key(occurrences$evidence_tier)),
      value = .plant_confidence_score(occurrences$confidence),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, data.frame(
      id = id,
      trait = paste0("occurrence_status__",
                     .plant_matrix_key(occurrences$occurrence_status)),
      value = suppressWarnings(as.numeric(occurrences$evidence_quality_score)),
      stringsAsFactors = FALSE
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "plant_part_group", "plant_part",
      exclude = c("unknown", "extract_unspecified")
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "tissue_group", "tissue",
      exclude = c("unknown", "extract_unspecified")
    ))
    long = rbind(long, .plant_context_matrix_long(
      occurrences, level, "method_group", "method",
      exclude = c("unknown")
    ))
  }
  trait_long = .plant_chemical_trait_matrix_long(
    occurrences = occurrences,
    categorate_result = categorate_result,
    level = level,
    profile = profile,
    min_confidence = min_confidence
  )
  if (nrow(trait_long) > 0) {
    long = rbind(long, trait_long)
  }
  long = long[!is.na(long$id) & long$id != "" &
                !is.na(long$trait) & long$trait != "compound__", ,
              drop = FALSE]
  .plant_wide_matrix(long, id_col = level, mode = mode,
                     max_traits = max_traits)
}

#' Filter plant phytochemistry occurrence evidence
#'
#' @description
#' Keeps plant-compound occurrence rows that meet explicit evidence, confidence,
#' and biological-context criteria. When `x` is a plant phytochemistry result,
#' the returned result is rebuilt with updated occurrence, compound-resolution,
#' trait-evidence, summary, matrix, and validation tables. Literature candidate
#' tables are preserved for audit/provenance.
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param occurrence_status Character vector of occurrence status values to
#' retain. Defaults to direct and curated reported records. Use `"all"` or
#' `NULL` to skip this filter.
#' @param analysis_ready Optional logical. If `TRUE`, retain only rows marked
#' analysis-ready; if `FALSE`, retain only non-analysis-ready rows; if `NULL`,
#' do not filter on the flag.
#' @param min_confidence Minimum confidence to retain. Use `NULL` to skip.
#' @param plant_part_group Optional plant-part groups to retain.
#' @param tissue_group Optional tissue groups to retain.
#' @param method_group Optional method groups to retain.
#' @param source_database Optional source databases to retain.
#' @param evidence_tier Optional evidence tiers to retain.
#'
#' @return A filtered plant phytochemistry result when `x` is a result object;
#' otherwise a filtered `PlantCompoundOccurrences` data frame.
#'
#' @export
filterPlantPhytochemistryEvidence = function(
    x,
    occurrence_status = c("direct_reported", "curated_reported"),
    analysis_ready = TRUE,
    min_confidence = "medium",
    plant_part_group = NULL,
    tissue_group = NULL,
    method_group = NULL,
    source_database = NULL,
    evidence_tier = NULL) {
  is_result = inherits(x, "uaf_plant_phytochemistry")
  occurrences = if (is_result) x$PlantCompoundOccurrences else x
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(if (is_result) .plant_rebuild_filtered_result(x, occurrences) else
      occurrences)
  }
  keep = rep(TRUE, nrow(occurrences))
  status_filter = .uaf_non_empty(occurrence_status)
  if (length(status_filter) > 0 && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$occurrence_status %in% status_filter
  }
  if (!is.null(analysis_ready)) {
    ready_value = .uaf_yes_no(isTRUE(analysis_ready))
    keep = keep & occurrences$analysis_ready == ready_value
  }
  if (!is.null(min_confidence)) {
    keep = keep & .plant_confidence_score(occurrences$confidence) >=
      .plant_confidence_score(min_confidence)
  }
  keep = .plant_apply_optional_filter(keep, occurrences$plant_part_group,
                                      plant_part_group)
  keep = .plant_apply_optional_filter(keep, occurrences$tissue_group,
                                      tissue_group)
  keep = .plant_apply_optional_filter(keep, occurrences$method_group,
                                      method_group)
  keep = .plant_apply_optional_filter(keep, occurrences$source_database,
                                      source_database)
  keep = .plant_apply_optional_filter(keep, occurrences$evidence_tier,
                                      evidence_tier)
  filtered = occurrences[keep, , drop = FALSE]
  row.names(filtered) = NULL
  if (is_result) .plant_rebuild_filtered_result(x, filtered) else filtered
}

#' Create a plant phytochemistry evidence review table
#'
#' @description
#' Builds a review worksheet for candidate, fallback, literature, or unresolved
#' plant-compound occurrence evidence. The table is intended for human review:
#' rows are not promoted automatically. Fill `review_decision` and optional
#' proposed fields, then pass the completed table to
#' `applyPlantPhytochemistryReview()`.
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param occurrence_status Occurrence statuses to include. Use `"all"` to
#' include every occurrence row.
#' @param include_analysis_ready Logical. If `FALSE`, rows already marked
#' analysis-ready are omitted unless `occurrence_status = "all"`.
#'
#' @return A data frame suitable for CSV export, editing, and re-import.
#'
#' @export
plantPhytochemistryReviewTable = function(
    x,
    occurrence_status = c("candidate", "taxon_fallback",
                          "literature_reported", "unresolved"),
    include_analysis_ready = FALSE) {
  occurrences = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$PlantCompoundOccurrences
  } else {
    x
  }
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) {
    return(.uaf_empty_table(.plant_review_cols()))
  }
  status_filter = .uaf_non_empty(occurrence_status)
  keep = rep(TRUE, nrow(occurrences))
  if (length(status_filter) > 0 && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$occurrence_status %in% status_filter
  }
  if (!isTRUE(include_analysis_ready) && !"all" %in% tolower(status_filter)) {
    keep = keep & occurrences$analysis_ready != "Yes"
  }
  review = occurrences[keep, , drop = FALSE]
  if (nrow(review) < 1) return(.uaf_empty_table(.plant_review_cols()))
  out = data.frame(
    review_id = sprintf("review_%04d", seq_len(nrow(review))),
    review_decision = "needs_review",
    reviewed_by = NA_character_,
    reviewed_at = NA_character_,
    review_note = NA_character_,
    proposed_evidence_tier = NA_character_,
    proposed_confidence = NA_character_,
    proposed_source_database = NA_character_,
    proposed_citation_or_url = NA_character_,
    proposed_plant_part = NA_character_,
    proposed_tissue = NA_character_,
    proposed_method = NA_character_,
    review_key = .plant_occurrence_key(review),
    review[, .plant_occurrence_cols(), drop = FALSE],
    stringsAsFactors = FALSE
  )
  row.names(out) = NULL
  out
}

#' Apply reviewed plant phytochemistry evidence decisions
#'
#' @description
#' Applies a completed review table created by
#' `plantPhytochemistryReviewTable()`. Supported decisions are
#' `"needs_review"`/`"keep"` (leave unchanged), `"keep_candidate"` (mark as
#' reviewed but keep candidate/fallback status), `"update_context"` (copy
#' proposed plant part/tissue/method/confidence/source fields),
#' `"promote_curated"` (promote to `manual_curated` evidence with review
#' provenance), and `"reject"`/`"exclude"` (remove the occurrence row).
#'
#' @param x Plant phytochemistry result or normalized occurrence table.
#' @param review_table Completed review table.
#' @param reviewer Optional reviewer name used when `reviewed_by` is empty.
#' @param require_citation Logical. If `TRUE`, `promote_curated` requires either
#' an existing evidence URL or `proposed_citation_or_url`.
#'
#' @return Updated plant phytochemistry result or occurrence table.
#'
#' @export
applyPlantPhytochemistryReview = function(x,
                                          review_table,
                                          reviewer = NA_character_,
                                          require_citation = TRUE) {
  if (!is.data.frame(review_table)) {
    stop("`review_table` must be a data frame.", call. = FALSE)
  }
  is_result = inherits(x, "uaf_plant_phytochemistry")
  occurrences = if (is_result) x$PlantCompoundOccurrences else x
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1 || nrow(review_table) < 1) {
    return(if (is_result) .plant_rebuild_filtered_result(x, occurrences) else
      occurrences)
  }
  names(review_table) = .plant_normalize_column_names(names(review_table))
  keep = rep(TRUE, nrow(occurrences))
  for (i in seq_len(nrow(review_table))) {
    review = review_table[i, , drop = FALSE]
    idx = .plant_review_match_indices(occurrences, review)
    if (length(idx) < 1) next
    decision = .plant_review_decision(review$review_decision)
    if (decision %in% c("reject", "exclude")) {
      keep[idx] = FALSE
      next
    }
    if (decision %in% c("needs_review", "keep", "")) next
    if (decision == "keep_candidate") {
      occurrences$curation_flag[idx] = "reviewed_candidate_kept"
      next
    }
    occurrences = .plant_apply_review_context(occurrences, idx, review)
    if (decision == "promote_curated") {
      citation = .uaf_first_non_empty_text(
        review$proposed_citation_or_url,
        occurrences$evidence_url[idx][[1]]
      )
      if (isTRUE(require_citation) && is.na(citation)) {
        warning("Skipping promotion for review row ", i,
                ": no citation or URL supplied.", call. = FALSE)
        next
      }
      occurrences$evidence_tier[idx] = "manual_curated"
      occurrences$confidence[idx] = .uaf_first_non_empty_text(
        review$proposed_confidence, "high"
      )
      occurrences$source_database[idx] = .uaf_first_non_empty_text(
        review$proposed_source_database, "manual_review"
      )
      occurrences$evidence_url[idx] = citation
      occurrences$curation_flag[idx] = "reviewed_curated"
      reviewer_value = .uaf_first_non_empty_text(review$reviewed_by,
                                                 reviewer)
      note = .uaf_first_non_empty_text(review$review_note)
      if (!is.na(reviewer_value) || !is.na(note)) {
        occurrences$evidence_text[idx] = .plant_truncate(
          .pubchem_collapse(c(occurrences$evidence_text[idx],
                              paste("reviewed_by", reviewer_value),
                              note)),
          1200
        )
      }
    }
  }
  updated = .plant_normalize_occurrences(occurrences[keep, , drop = FALSE])
  updated = .plant_collapse_occurrence_evidence(updated)
  if (!is_result) return(updated)
  out = .plant_rebuild_filtered_result(x, updated)
  out$Provenance = .plant_bind_tables(list(
    out$Provenance,
    .plant_provenance("evidence_review", "manual_review",
                      paste(unique(.uaf_non_empty(updated$species)),
                            collapse = "; "),
                      NA_character_, nrow(review_table),
                      "Human review table applied to plant occurrence evidence.")
  ), .plant_provenance_cols())
  out$Validation = validatePlantPhytochemistryResult(out)
  out
}

#' Join plant chemistry summaries to metadata
#'
#' @param metadata User metadata data frame.
#' @param plant_chemistry Plant phytochemistry result, summary table, or matrix.
#' @param species_col Species column in `metadata`.
#' @param join_level One of `"species"`, `"genus"`, or `"family"`.
#'
#' @return Joined data frame.
#'
#' @export
joinPlantChemistryMetadata = function(metadata,
                                      plant_chemistry,
                                      species_col = "species",
                                      join_level = "species") {
  if (!is.data.frame(metadata)) stop("`metadata` must be a data frame.",
                                     call. = FALSE)
  if (!species_col %in% names(metadata)) {
    stop("`species_col` was not found in metadata.", call. = FALSE)
  }
  chemistry = if (inherits(plant_chemistry, "uaf_plant_phytochemistry")) {
    plant_chemistry$SpeciesChemistrySummary
  } else {
    plant_chemistry
  }
  if (!is.data.frame(chemistry)) {
    stop("`plant_chemistry` must be a result object or data frame.",
         call. = FALSE)
  }
  if (!join_level %in% names(chemistry)) {
    stop("`join_level` was not found in the chemistry table.", call. = FALSE)
  }
  metadata$.uaf_join_key = .plant_clean_name(metadata[[species_col]])
  chemistry$.uaf_join_key = .plant_clean_name(chemistry[[join_level]])
  out = merge(metadata, chemistry, by = ".uaf_join_key", all.x = TRUE)
  out$.uaf_join_key = NULL
  row.names(out) = NULL
  out
}

#' Score plant chemistry candidates
#'
#' @param x Plant phytochemistry result or species summary table.
#' @param weights Named numeric vector for transparent score components.
#'
#' @return Data frame with component scores and total prioritization score.
#'
#' @export
scorePlantChemistryCandidates = function(x,
                                         weights = c(direct_species = 0.25,
                                                     resolved_compounds = 0.20,
                                                     source_coverage = 0.20,
                                                     literature = 0.15,
                                                     traits = 0.20)) {
  summary = if (inherits(x, "uaf_plant_phytochemistry")) {
    x$SpeciesChemistrySummary
  } else {
    x
  }
  if (!is.data.frame(summary) || nrow(summary) < 1) {
    return(.uaf_empty_table(c("species", "chemistry_priority_score")))
  }
  weights = weights[!is.na(weights)]
  if (length(weights) < 1 || sum(weights) <= 0) {
    weights = c(direct_species = 1)
  }
  weights = weights / sum(weights)
  norm = function(value) {
    value = suppressWarnings(as.numeric(value))
    max_value = max(value, na.rm = TRUE)
    if (!is.finite(max_value) || max_value <= 0) return(rep(0, length(value)))
    pmin(value / max_value, 1)
  }
  components = data.frame(
    species = summary$species,
    direct_species_score = norm(summary$direct_species_compound_count),
    resolved_compound_score = norm(summary$resolved_compound_count),
    source_coverage_score = suppressWarnings(as.numeric(summary$source_coverage_score)),
    literature_score = norm(summary$literature_candidate_count),
    trait_score = norm(summary$chemical_trait_count),
    stringsAsFactors = FALSE
  )
  components$source_coverage_score[is.na(components$source_coverage_score)] = 0
  total = rep(0, nrow(components))
  if ("direct_species" %in% names(weights)) {
    total = total + weights[["direct_species"]] * components$direct_species_score
  }
  if ("resolved_compounds" %in% names(weights)) {
    total = total + weights[["resolved_compounds"]] * components$resolved_compound_score
  }
  if ("source_coverage" %in% names(weights)) {
    total = total + weights[["source_coverage"]] * components$source_coverage_score
  }
  if ("literature" %in% names(weights)) {
    total = total + weights[["literature"]] * components$literature_score
  }
  if ("traits" %in% names(weights)) {
    total = total + weights[["traits"]] * components$trait_score
  }
  components$chemistry_priority_score = round(total, 4)
  components$score_interpretation = "Prioritization aid only; not evidence of efficacy, occurrence completeness, remediation performance, or causal mechanism."
  components
}

#' Export plant phytochemistry result tables
#'
#' @param x Plant phytochemistry result.
#' @param path Output directory for `format = "csv"` or `.xlsx` path for
#' `format = "xlsx"`.
#' @param format One of `"csv"`, `"xlsx"`, or `"auto"`.
#' @param tables Optional table names to export.
#' @param preset One of `"all"` or `"analysis_ready"`. The latter filters to
#' direct/curated analysis-ready occurrence evidence before export while
#' preserving validation and provenance tables.
#' @param include_empty Logical. If `TRUE`, include empty schema tables.
#' @param overwrite Logical. If `TRUE`, replace existing output.
#' @param max_cell_chars Maximum characters retained in any character cell.
#'
#' @return Export manifest data frame.
#'
#' @export
exportPlantPhytochemistryWorkbook = function(x,
                                             path,
                                             format = c("auto", "xlsx", "csv"),
                                             tables = NULL,
                                             preset = c("all", "analysis_ready"),
                                             include_empty = TRUE,
                                             overwrite = FALSE,
                                             max_cell_chars = 30000) {
  format = match.arg(format)
  preset = match.arg(preset)
  if (missing(path) || is.null(path) || length(.uaf_non_empty(path)) < 1) {
    stop("`path` is required.", call. = FALSE)
  }
  if (!is.list(x) || is.data.frame(x)) {
    stop("`x` must be a plant phytochemistry result list.", call. = FALSE)
  }
  if (preset == "analysis_ready") {
    x = filterPlantPhytochemistryEvidence(
      x,
      occurrence_status = c("direct_reported", "curated_reported"),
      analysis_ready = TRUE,
      min_confidence = "medium"
    )
  }
  resolved = .categorate_export_resolve_path(path, format)
  export_tables = .plant_export_tables(x, tables, include_empty, max_cell_chars)
  if (length(export_tables) < 1) stop("No exportable tables found.",
                                      call. = FALSE)

  all_names = c("ExportManifest", names(export_tables))
  sheet_names = .categorate_export_sheet_names(all_names)
  file_names = .categorate_export_file_names(all_names)
  created_at = .plant_timestamp()
  manifest = data.frame(
    Table = all_names,
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
  }
  row.names(manifest) = NULL
  manifest
}

#' @export
print.uaf_plant_phytochemistry = function(x, ...) {
  cat("uafR plant phytochemistry result\n")
  cat("  plant queries: ", nrow(x$PlantQueries), "\n", sep = "")
  cat("  occurrence rows: ", nrow(x$PlantCompoundOccurrences), "\n", sep = "")
  cat("  literature candidates: ", nrow(x$LiteratureCandidates), "\n", sep = "")
  cat("  resolved compounds: ",
      sum(x$CompoundResolution$resolved %in% TRUE), "\n", sep = "")
  if (is.list(x$Validation) && is.data.frame(x$Validation$Summary)) {
    cat("  validation status: ", x$Validation$Summary$Status[[1]], "\n",
        sep = "")
  }
  invisible(x)
}

.plant_data_dictionary = function() {
  specs = list(
    .plant_schema("PlantQueries", .plant_query_cols(),
                  required = c("query_id", "query_plant",
                               "query_plant_clean", "species", "species_slug"),
                  role = "query",
                  description = "Input plant names and query-normalized fields."),
    .plant_schema("PlantNameResolution", .plant_name_resolution_cols(),
                  required = c("query_id", "query_plant",
                               "query_status", "matched_taxon"),
                  role = "taxonomy",
                  description = "Plant name parsing and resolution status."),
    .plant_schema("ProviderDiagnostics", .plant_provider_diagnostic_cols(),
                  required = c("provider", "enabled", "queried",
                               "record_count", "error_count", "retrieved_at"),
                  role = "diagnostics",
                  description = "Provider availability, requests, records, and warnings."),
    .plant_schema("PlantCompoundOccurrences", .plant_occurrence_cols(),
                  required = c("query_plant", "query_plant_clean",
                               "species", "compound_name",
                               "compound_name_clean", "source_database",
                               "confidence", "evidence_tier",
                               "occurrence_status", "occurrence_basis",
                               "analysis_ready", "evidence_quality_score"),
                  role = "occurrence",
                  allowed = list(
                    confidence = .plant_confidence_values(),
                    evidence_tier = .plant_evidence_tiers(),
                    plant_part_group = .plant_context_group_values(),
                    tissue_group = .plant_context_group_values(),
                    method_group = .plant_method_group_values(),
                    biological_context_status = .plant_context_status_values(),
                    occurrence_status = .plant_occurrence_status_values(),
                    occurrence_basis = .plant_occurrence_basis_values(),
                    analysis_ready = .plant_yes_no_values(),
                    matched_rank = c("species", "genus", "family",
                                     "unknown")
                  ),
                  description = "Normalized species-compound occurrence evidence."),
    .plant_schema("LiteratureCandidates", .plant_literature_cols(),
                  required = c("query_plant", "species", "source_database",
                               "evidence_tier", "confidence"),
                  role = "literature",
                  allowed = list(confidence = .plant_confidence_values(),
                                 evidence_tier = .plant_evidence_tiers()),
                  description = "Publication and annotation candidates; not confirmed occurrence by default."),
    .plant_schema("CompoundResolution", .plant_compound_resolution_cols(),
                  required = c("compound_name", "compound_name_clean",
                               "resolved"),
                  role = "compound_identity",
                  description = "Compound identity resolution from uafR enrichment."),
    .plant_schema("SpeciesChemistrySummary", .plant_summary_cols(),
                  required = c("species", "species_slug", "compound_count",
                               "source_coverage_score"),
                  role = "summary",
                  description = "One row per species with counts, coverage, traits, and caveats."),
    .plant_schema("SpeciesChemistryMatrix", c("species"),
                  required = c("species"),
                  role = "matrix",
                  description = "Wide species-feature matrix; columns vary by profile and data."),
    .plant_schema("TraitEvidence", .plant_trait_evidence_cols(),
                  required = c("species", "compound_name_clean",
                               "EvidenceType", "AnalysisKey",
                               "SourceDatabase", "Confidence"),
                  role = "evidence",
                  description = "Species-compound links to uafR chemical trait evidence."),
    .plant_schema("Validation", c("Summary", "TableQuality",
                                  "ProviderDiagnostics", "Issues",
                                  "DataDictionary"),
                  required = c("Summary", "TableQuality", "Issues"),
                  role = "validation",
                  description = "Nested validation result."),
    .plant_schema("DataDictionary", c("Table", "Column", "ExpectedType",
                                      "Required", "Role", "AllowedValues",
                                      "Description"),
                  required = c("Table", "Column", "Required"),
                  role = "dictionary",
                  description = "Schema table."),
    .plant_schema("Provenance", .plant_provenance_cols(),
                  required = c("step", "source", "retrieved_at",
                               "record_count"),
                  role = "provenance",
                  description = "Source, retrieval, and assembly provenance.")
  )
  out = do.call(rbind, specs)
  row.names(out) = NULL
  out
}

.plant_schema = function(table, cols, required = character(), role,
                         description, allowed = list()) {
  data.frame(
    Table = table,
    Column = cols,
    ExpectedType = .plant_expected_type(cols),
    Required = cols %in% required,
    Role = role,
    AllowedValues = vapply(cols, function(col) {
      if (!is.null(allowed[[col]])) paste(allowed[[col]], collapse = "; ") else ""
    }, character(1)),
    Description = description,
    stringsAsFactors = FALSE
  )
}

.plant_expected_type = function(cols) {
  integer_cols = grepl(
    paste(c("^input_order$", "_count$", "^RowCount$", "^ColumnCount$",
            "^request_count$", "^cache_hit_count$", "^record_count$",
            "^error_count$", "^warning_count$", "^max_pubmed_records$",
            "^max_provider_records$"),
          collapse = "|"),
    cols,
    ignore.case = TRUE
  )
  ifelse(integer_cols,
         "integer",
         ifelse(grepl("score|fraction", cols, ignore.case = TRUE),
                "numeric", "character"))
}

.plant_query_cols = function() {
  c("query_id", "input_order", "query_plant", "query_plant_clean",
    "species", "genus", "family", "species_slug", "taxon_fallback")
}

.plant_name_resolution_cols = function() {
  c("query_id", "query_plant", "query_plant_clean", "matched_taxon",
    "matched_rank", "species", "genus", "family", "species_slug",
    "query_status", "resolution_source", "resolution_note")
}

.plant_provider_diagnostic_cols = function() {
  c("provider", "enabled", "queried", "available", "request_count",
    "cache_hit_count", "record_count", "error_count", "warning_count",
    "status", "message", "retrieved_at")
}

.plant_occurrence_cols = function() {
  c("query_plant", "query_plant_clean", "matched_taxon", "matched_rank",
    "species", "genus", "family", "compound_name", "compound_name_clean",
    "compound_id", "compound_id_type", "source_database",
    "source_record_id", "evidence_text", "evidence_url", "reference_id",
    "pmid", "doi", "plant_part", "plant_part_group", "tissue",
    "tissue_group", "method", "method_group", "biological_context_status",
    "occurrence_type", "occurrence_status", "occurrence_basis",
    "analysis_ready", "evidence_quality_score", "retrieved_at",
    "confidence", "curation_flag", "evidence_tier")
}

.plant_literature_cols = function() {
  c("query_plant", "query_plant_clean", "species", "genus", "family",
    "source_database", "source_record_id", "pmid", "doi", "title",
    "abstract", "chemical_mention", "species_mention", "evidence_text",
    "evidence_url", "retrieved_at", "confidence", "curation_flag",
    "evidence_tier")
}

.plant_compound_resolution_cols = function() {
  c("compound_name", "compound_name_clean", "query_count", "resolved", "CID",
    "InChIKey", "SMILES", "MolecularFormula", "resolution_source", "notes")
}

.plant_summary_cols = function() {
  c("species", "species_slug", "genus", "family", "query_status",
    "compound_count", "resolved_compound_count", "unresolved_compound_count",
    "analysis_ready_compound_count",
    "direct_species_compound_count", "genus_level_compound_count",
    "family_level_compound_count", "literature_candidate_count",
    "database_occurrence_count", "direct_reported_occurrence_count",
    "curated_reported_occurrence_count", "candidate_occurrence_count",
    "taxon_fallback_occurrence_count", "context_known_occurrence_count",
    "mean_evidence_quality_score", "analysis_ready_fraction",
    "plant_part_groups", "tissue_groups", "method_groups",
    "evidence_tier_summary",
    "source_database_count", "source_databases",
    "natural_product_superclasses", "natural_product_classes",
    "natural_product_subclasses", "dominant_compound_classes",
    "kingdoms_observed", "families_observed", "is_plant_occurring",
    "volatile_proxy_fraction", "lipophilic_fraction", "oxygenated_fraction",
    "nitrogenous_fraction", "sulfur_containing_fraction",
    "halogenated_fraction", "kegg_pathway_groups", "chemical_trait_count",
    "high_confidence_trait_count", "source_coverage_score",
    "uafR_validation_status")
}

.plant_trait_evidence_cols = function() {
  c("species", "species_slug", "compound_name", "compound_name_clean",
    "EvidenceType", "AnalysisKey", "SourceDatabase", "EvidenceText",
    "EvidenceURL", "Confidence")
}

.plant_provenance_cols = function() {
  c("step", "source", "query", "source_url", "cache_file", "retrieved_at",
    "record_count", "notes")
}

.plant_validation_issue_cols = function() {
  c("severity", "table", "column", "issue", "expected", "observed",
    "row_count", "examples")
}

.plant_table_quality_cols = function() {
  c("Table", "Present", "RowCount", "ColumnCount", "MissingRequiredColumns",
    "MissingOptionalColumns", "DuplicateKeyColumns", "Status")
}

.plant_confidence_values = function() c("low", "medium", "high", "unknown")

.plant_yes_no_values = function() c("Yes", "No")

.plant_context_group_values = function() {
  c("unknown", "root_belowground", "leaf", "stem_shoot", "bark_wood",
    "flower", "fruit_seed", "aerial", "whole_plant",
    "exudate_rhizosphere", "vascular", "secretory", "epidermal",
    "microbial", "extract_unspecified", "other")
}

.plant_method_group_values = function() {
  c("unknown", "gc_ms", "lc_ms", "hplc", "nmr", "mass_spectrometry",
    "chromatography", "spectroscopy", "database_record",
    "literature_curation", "other")
}

.plant_occurrence_status_values = function() {
  c("direct_reported", "curated_reported", "literature_reported",
    "candidate", "taxon_fallback", "unresolved", "unknown")
}

.plant_occurrence_basis_values = function() {
  c("species_database_record", "species_curated_record",
    "species_literature_record", "candidate_co_mention",
    "genus_fallback_record", "family_fallback_record", "unresolved",
    "unknown")
}

.plant_context_status_values = function() {
  c("plant_part_and_method_known", "plant_part_known", "method_known",
    "context_missing")
}

.plant_evidence_tiers = function() {
  c("direct_species_database", "direct_species_literature",
    "direct_species_pubtator_candidate", "genus_database_fallback",
    "genus_literature_fallback", "family_database_fallback",
    "manual_curated", "unresolved")
}

.plant_queries = function(plants, taxon_fallback = c("species", "genus")) {
  if (missing(plants) || is.null(plants)) {
    stop("`plants` must contain at least one plant name.", call. = FALSE)
  }
  if (is.data.frame(plants)) {
    names(plants) = .plant_normalize_column_names(names(plants))
    species = if ("species" %in% names(plants)) plants$species else plants[[1]]
    genus = if ("genus" %in% names(plants)) plants$genus else .plant_genus(species)
    family = if ("family" %in% names(plants)) plants$family else
      rep(NA_character_, length(species))
  } else {
    species = plants
    genus = .plant_genus(species)
    family = rep(NA_character_, length(species))
  }
  species = .uaf_squish_text(species)
  keep = !is.na(species) & species != ""
  species = species[keep]
  genus = .uaf_squish_text(genus)[keep]
  family = .uaf_squish_text(family)[keep]
  if (length(species) < 1) {
    stop("`plants` must contain at least one non-empty plant name.",
         call. = FALSE)
  }
  data.frame(
    query_id = paste0("plant_", seq_along(species)),
    input_order = seq_along(species),
    query_plant = species,
    query_plant_clean = .plant_clean_name(species),
    species = species,
    genus = genus,
    family = family,
    species_slug = .plant_slug(species),
    taxon_fallback = paste(.uaf_non_empty(taxon_fallback), collapse = "; "),
    stringsAsFactors = FALSE
  )
}

.plant_name_resolution = function(plant_queries) {
  out = data.frame(
    query_id = plant_queries$query_id,
    query_plant = plant_queries$query_plant,
    query_plant_clean = plant_queries$query_plant_clean,
    matched_taxon = plant_queries$species,
    matched_rank = "species",
    species = plant_queries$species,
    genus = plant_queries$genus,
    family = plant_queries$family,
    species_slug = plant_queries$species_slug,
    query_status = ifelse(.plant_is_binomial(plant_queries$species),
                          "parsed_species", "name_needs_review"),
    resolution_source = "input_name_parser",
    resolution_note = ifelse(.plant_is_binomial(plant_queries$species),
                             "Species-like binomial parsed from input.",
                             "Input is not a clear binomial; review taxon name."),
    stringsAsFactors = FALSE
  )
  row.names(out) = NULL
  out
}

.plant_normalize_sources = function(sources) {
  allowed = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator")
  if (is.null(sources)) return(allowed)
  sources = tolower(.uaf_non_empty(sources))
  if (length(sources) < 1) return(character())
  bad = setdiff(sources, allowed)
  if (length(bad) > 0) {
    warning("Unsupported plant provider(s) skipped: ",
            paste(bad, collapse = ", "), call. = FALSE)
  }
  unique(intersect(sources, allowed))
}

.plant_provider_dispatch = function(plant_queries, sources, cache, cache_dir,
                                    throttle, ncbi_email, ncbi_tool,
                                    ncbi_api_key, max_pubmed_records,
                                    max_provider_records, provider_results,
                                    request_fun, pubtator_request_fun,
                                    refresh) {
  occurrence_rows = list()
  literature_rows = list()
  diagnostic_rows = list()
  provenance_rows = list()
  for (provider in sources) {
    result = .plant_query_provider(
      provider = provider,
      plant_queries = plant_queries,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      ncbi_email = ncbi_email,
      ncbi_tool = ncbi_tool,
      ncbi_api_key = ncbi_api_key,
      max_pubmed_records = max_pubmed_records,
      max_provider_records = max_provider_records,
      provider_results = provider_results,
      request_fun = request_fun,
      pubtator_request_fun = pubtator_request_fun,
      refresh = refresh
    )
    occurrence_rows[[length(occurrence_rows) + 1]] =
      result$PlantCompoundOccurrences
    literature_rows[[length(literature_rows) + 1]] =
      result$LiteratureCandidates
    diagnostic_rows[[length(diagnostic_rows) + 1]] =
      result$ProviderDiagnostics
    provenance_rows[[length(provenance_rows) + 1]] = result$Provenance
  }
  list(
    PlantCompoundOccurrences = .plant_bind_occurrences(occurrence_rows),
    LiteratureCandidates = .plant_bind_tables(literature_rows,
                                              .plant_literature_cols()),
    ProviderDiagnostics = .plant_bind_tables(diagnostic_rows,
                                             .plant_provider_diagnostic_cols()),
    Provenance = .plant_bind_tables(provenance_rows, .plant_provenance_cols())
  )
}

.plant_query_provider = function(provider, plant_queries, cache, cache_dir,
                                 throttle, ncbi_email, ncbi_tool,
                                 ncbi_api_key, max_pubmed_records,
                                 max_provider_records, provider_results,
                                 request_fun, pubtator_request_fun, refresh) {
  if (!is.null(provider_results) && !is.null(provider_results[[provider]])) {
    return(.plant_provider_from_result(provider, provider_results[[provider]],
                                       plant_queries))
  }
  switch(provider,
         lotus = .plant_query_lotus(
           plant_queries, cache, cache_dir, throttle, request_fun,
           max_provider_records),
         knapsack = .plant_query_knapsack(
           plant_queries, cache, cache_dir, throttle, request_fun,
           max_provider_records),
         npass = .plant_query_npass(
           plant_queries),
         pubchem = .plant_query_pubchem_occurrences(
           plant_queries, cache, cache_dir, throttle, ncbi_email, ncbi_tool,
           ncbi_api_key, max_provider_records, request_fun),
         pubmed = .plant_query_pubmed_literature(
           plant_queries, cache, cache_dir, throttle, ncbi_email, ncbi_tool,
           ncbi_api_key, max_pubmed_records, request_fun),
         pubtator = .plant_query_pubtator_literature(
           plant_queries, cache, cache_dir, throttle, pubtator_request_fun),
         .plant_provider_empty_result(provider, plant_queries,
                                      status = "not_implemented",
                                      message = paste0(
                                        provider,
                                        " species-first adapter is scaffolded ",
                                        "but requires a public download/API ",
                                        "parser or provider_results input."))
  )
}

.plant_provider_from_result = function(provider, result, plant_queries) {
  if (is.data.frame(result)) {
    occurrences = .plant_normalize_occurrences(result, source_hint = provider)
    literature = .uaf_empty_table(.plant_literature_cols())
  } else if (is.list(result)) {
    occurrences = .plant_normalize_occurrences(
      result$PlantCompoundOccurrences %||% result$occurrences %||%
        .uaf_empty_table(.plant_occurrence_cols()),
      source_hint = provider
    )
    literature = .plant_normalize_literature(
      result$LiteratureCandidates %||% result$literature_candidates %||%
        .uaf_empty_table(.plant_literature_cols()),
      source_hint = provider
    )
  } else {
    occurrences = .plant_empty_occurrences()
    literature = .uaf_empty_table(.plant_literature_cols())
  }
  occurrences = .plant_match_occurrences_to_queries(occurrences, plant_queries)
  literature = .plant_match_literature_to_queries(literature, plant_queries)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = literature,
    ProviderDiagnostics = .plant_provider_diagnostics(
      provider, TRUE, TRUE, TRUE, 0, 0,
      nrow(occurrences) + nrow(literature), 0, 0, "ok",
      "Provider rows supplied by provider_results."),
    Provenance = .plant_provenance("provider_dispatch", provider,
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   NA_character_,
                                   nrow(occurrences) + nrow(literature),
                                   "Provider rows supplied locally.")
  )
}

.plant_query_lotus = function(plant_queries, cache, cache_dir, throttle,
                              request_fun, max_records) {
  rows = list()
  request_count = 0
  errors = 0
  for (i in seq_len(nrow(plant_queries))) {
    plant = plant_queries$species[[i]]
    url = paste0("https://lotus.naturalproducts.net/api/search/simple?query=",
                 utils::URLencode(plant, reserved = TRUE))
    request_count = request_count + 1
    result = tryCatch({
      .plant_fetch_json(url, cache, file.path(cache_dir, "lotus"),
                        throttle, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    parsed = .plant_lotus_rows(plant_queries[i, , drop = FALSE], result,
                               url, max_records)
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "lotus", TRUE, TRUE, TRUE, request_count, NA_integer_,
      nrow(occurrences), errors, 0, status,
      "LOTUS simple API taxon-oriented candidates; only rows with taxon evidence are retained."),
    Provenance = .plant_provenance("provider_dispatch", "lotus",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://lotus.naturalproducts.net/api/search/simple",
                                   nrow(occurrences),
                                   "LOTUS Natural Products Online candidate search.")
  )
}

.plant_query_knapsack = function(plant_queries, cache, cache_dir, throttle,
                                 request_fun, max_records) {
  rows = list()
  request_count = 0
  errors = 0
  for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    terms = .plant_taxon_query_terms(query_row)
    for (j in seq_len(nrow(terms))) {
      url = paste0(
        "https://www.knapsackfamily.com/knapsack_core/result.php?sname=organism&word=",
        utils::URLencode(terms$term[[j]], reserved = TRUE)
      )
      request_count = request_count + 1
      html = tryCatch({
        .plant_fetch_text(url, cache, file.path(cache_dir, "knapsack"),
                          throttle, request_fun)
      }, error = function(error) {
        errors <<- errors + 1
        NULL
      })
      parsed = .plant_knapsack_rows(query_row, html, terms$term[[j]],
                                    terms$rank[[j]], url, max_records)
      if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
    }
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "knapsack", TRUE, TRUE, TRUE, request_count, NA_integer_,
      nrow(occurrences), errors, 0, status,
      "KNApSAcK organism query records normalized from public HTML output."),
    Provenance = .plant_provenance("provider_dispatch", "knapsack",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://www.knapsackfamily.com/knapsack_core/result.php",
                                   nrow(occurrences),
                                   "KNApSAcK organism-metabolite candidate search.")
  )
}

.plant_query_npass = function(plant_queries) {
  .plant_provider_empty_result(
    "npass", plant_queries, status = "not_queried",
    message = paste("NPASS species-source data are exposed primarily through",
                    "web search and downloadable files. Supply NPASS rows via",
                    "`provider_results` or curated intake until a small stable",
                    "species-query endpoint is added."))
}

.plant_query_pubchem_occurrences = function(plant_queries, cache, cache_dir,
                                            throttle, ncbi_email, ncbi_tool,
                                            ncbi_api_key, max_records,
                                            request_fun) {
  rows = list()
  request_count = 0
  errors = 0
  effective_throttle = throttle
  if (is.null(request_fun)) {
    has_ncbi_key = length(.uaf_non_empty(ncbi_api_key)) > 0
    effective_throttle = max(throttle, ifelse(has_ncbi_key, 0.10, 0.34))
  }
  for (i in seq_len(nrow(plant_queries))) {
    query_row = plant_queries[i, , drop = FALSE]
    request_count = request_count + 1
    taxid_search = tryCatch({
      .plant_ncbi_taxonomy_search(query_row$species[[1]], cache, cache_dir,
                                  effective_throttle, ncbi_email, ncbi_tool,
                                  ncbi_api_key, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    taxids = .plant_pubmed_ids(taxid_search)
    if (length(taxids) < 1) next
    taxid = taxids[[1]]
    url = paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/taxonomy/",
                 utils::URLencode(taxid, reserved = TRUE), "/JSON")
    request_count = request_count + 1
    result = tryCatch({
      .plant_fetch_json(url, cache, file.path(cache_dir, "pubchem_taxonomy"),
                        effective_throttle, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    parsed = .plant_pubchem_taxonomy_rows(query_row, taxid, result, url,
                                          max_records)
    if (nrow(parsed) > 0) rows[[length(rows) + 1]] = parsed
  }
  occurrences = .plant_bind_occurrences(rows)
  status = .plant_provider_status(nrow(occurrences), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubchem", TRUE, TRUE, TRUE, request_count, NA_integer_,
      nrow(occurrences), errors, 0, status,
      "PubChem taxonomy PUG-View chemical annotations; review evidence context before treating as occurrence."),
    Provenance = .plant_provenance("provider_dispatch", "pubchem",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/taxonomy/",
                                   nrow(occurrences),
                                   "PubChem taxonomy annotation search.")
  )
}

.plant_query_pubmed_literature = function(plant_queries, cache, cache_dir,
                                          throttle, ncbi_email, ncbi_tool,
                                          ncbi_api_key, max_pubmed_records,
                                          request_fun) {
  rows = list()
  request_count = 0
  errors = 0
  effective_throttle = throttle
  if (is.null(request_fun)) {
    has_ncbi_key = length(.uaf_non_empty(ncbi_api_key)) > 0
    effective_throttle = max(throttle, ifelse(has_ncbi_key, 0.10, 0.34))
  }
  for (i in seq_len(nrow(plant_queries))) {
    plant = plant_queries$query_plant[[i]]
    term = .plant_pubmed_query(plant)
    esearch_url = .plant_ncbi_url(
      endpoint = "esearch.fcgi",
      params = c(db = "pubmed", term = term, retmode = "json",
                 retmax = as.character(max_pubmed_records),
                 tool = ncbi_tool, email = ncbi_email),
      api_key = ncbi_api_key
    )
    request_count = request_count + 1
    search = tryCatch({
      .plant_fetch_json(esearch_url, cache, file.path(cache_dir, "pubmed"),
                        effective_throttle, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    ids = .plant_pubmed_ids(search)
    if (length(ids) < 1) next
    summary_url = .plant_ncbi_url(
      endpoint = "esummary.fcgi",
      params = c(db = "pubmed", id = paste(ids, collapse = ","),
                 retmode = "json", tool = ncbi_tool, email = ncbi_email),
      api_key = ncbi_api_key
    )
    request_count = request_count + 1
    summary = tryCatch({
      .plant_fetch_json(summary_url, cache, file.path(cache_dir, "pubmed"),
                        effective_throttle, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    rows[[length(rows) + 1]] =
      .plant_pubmed_summary_rows(plant_queries[i, , drop = FALSE], ids,
                                 summary)
  }
  literature = .plant_bind_tables(rows, .plant_literature_cols())
  status = .plant_provider_status(nrow(literature), errors)
  list(
    PlantCompoundOccurrences = .plant_empty_occurrences(),
    LiteratureCandidates = literature,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubmed", TRUE, TRUE, TRUE, request_count, NA_integer_,
      nrow(literature), errors, 0, status,
      "PubMed E-utilities literature candidates; not confirmed occurrence."),
    Provenance = .plant_provenance("provider_dispatch", "pubmed",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                                   nrow(literature),
                                   "PubMed candidate literature search.")
  )
}

.plant_query_pubtator_literature = function(plant_queries, cache, cache_dir,
                                            throttle, request_fun) {
  rows = list()
  request_count = 0
  errors = 0
  effective_throttle = if (is.null(request_fun)) max(throttle, 0.34) else
    throttle
  for (i in seq_len(nrow(plant_queries))) {
    plant = plant_queries$query_plant[[i]]
    url = paste0("https://www.ncbi.nlm.nih.gov/research/pubtator3-api/search/",
                 "?text=", utils::URLencode(
                   paste(plant, "phytochemical metabolite"),
                   reserved = TRUE))
    request_count = request_count + 1
    result = tryCatch({
      .plant_fetch_json(url, cache, file.path(cache_dir, "pubtator"),
                        effective_throttle, request_fun)
    }, error = function(error) {
      errors <<- errors + 1
      NULL
    })
    parsed = .plant_pubtator_rows(plant_queries[i, , drop = FALSE], result)
    rows[[length(rows) + 1]] = parsed$LiteratureCandidates
  }
  literature = .plant_bind_tables(rows, .plant_literature_cols())
  occurrences = .plant_occurrences_from_pubtator(literature)
  status = .plant_provider_status(nrow(literature), errors)
  list(
    PlantCompoundOccurrences = occurrences,
    LiteratureCandidates = literature,
    ProviderDiagnostics = .plant_provider_diagnostics(
      "pubtator", TRUE, TRUE, TRUE, request_count, NA_integer_,
      nrow(literature), errors, 0, status,
      "PubTator candidate co-mentions; not confirmed occurrence."),
    Provenance = .plant_provenance("provider_dispatch", "pubtator",
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   "https://www.ncbi.nlm.nih.gov/research/pubtator3-api/",
                                   nrow(literature),
                                   "PubTator chemical/species candidate search.")
  )
}

.plant_provider_empty_result = function(provider, plant_queries,
                                        status = "no_records",
                                        message = "No provider records returned.") {
  list(
    PlantCompoundOccurrences = .plant_empty_occurrences(),
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    ProviderDiagnostics = .plant_provider_diagnostics(
      provider = provider,
      enabled = TRUE,
      queried = FALSE,
      available = FALSE,
      request_count = 0,
      cache_hit_count = 0,
      record_count = 0,
      error_count = 0,
      warning_count = ifelse(status %in% c("not_implemented", "not_queried"),
                             1, 0),
      status = status,
      message = message
    ),
    Provenance = .plant_provenance("provider_dispatch", provider,
                                   paste(plant_queries$query_plant,
                                         collapse = "; "),
                                   NA_character_, 0, message)
  )
}

.plant_provider_status = function(record_count, errors) {
  if (errors > 0) return("warning")
  if (record_count > 0) return("ok")
  "no_records"
}

.plant_provider_diagnostics = function(provider, enabled, queried, available,
                                       request_count, cache_hit_count,
                                       record_count, error_count,
                                       warning_count, status, message) {
  data.frame(
    provider = provider,
    enabled = .uaf_yes_no(enabled),
    queried = .uaf_yes_no(queried),
    available = .uaf_yes_no(available),
    request_count = suppressWarnings(as.integer(request_count)),
    cache_hit_count = suppressWarnings(as.integer(cache_hit_count)),
    record_count = suppressWarnings(as.integer(record_count)),
    error_count = suppressWarnings(as.integer(error_count)),
    warning_count = suppressWarnings(as.integer(warning_count)),
    status = status,
    message = message,
    retrieved_at = .plant_timestamp(),
    stringsAsFactors = FALSE
  )
}

.plant_normalize_occurrences = function(x, source_hint = NA_character_) {
  cols = .plant_occurrence_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_table(cols))
  names(x) = .plant_normalize_column_names(names(x))
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  x$query_plant = .uaf_squish_text(x$query_plant)
  x$species = .uaf_squish_text(x$species)
  missing_query = is.na(x$query_plant) | x$query_plant == ""
  x$query_plant[missing_query] = x$species[missing_query]
  x$query_plant_clean = .plant_clean_name(
    ifelse(is.na(x$query_plant_clean) | x$query_plant_clean == "",
           x$query_plant, x$query_plant_clean)
  )
  x$genus = ifelse(is.na(x$genus) | x$genus == "", .plant_genus(x$species),
                   .uaf_squish_text(x$genus))
  x$family = .uaf_squish_text(x$family)
  x$matched_taxon = .uaf_squish_text(ifelse(is.na(x$matched_taxon) |
                                              x$matched_taxon == "",
                                            x$species, x$matched_taxon))
  x$matched_rank = tolower(.uaf_squish_text(ifelse(is.na(x$matched_rank) |
                                                     x$matched_rank == "",
                                                   "species", x$matched_rank)))
  x$matched_rank[!x$matched_rank %in% c("species", "genus", "family")] = "unknown"
  x$compound_name = .uaf_squish_text(x$compound_name)
  x$compound_name_clean = .plant_clean_compound(
    ifelse(is.na(x$compound_name_clean) | x$compound_name_clean == "",
           x$compound_name, x$compound_name_clean)
  )
  x$source_database = .uaf_squish_text(ifelse(is.na(x$source_database) |
                                                x$source_database == "",
                                              source_hint, x$source_database))
  x$source_database[is.na(x$source_database) | x$source_database == ""] =
    "unknown"
  x$retrieved_at = ifelse(is.na(x$retrieved_at) | x$retrieved_at == "",
                          .plant_timestamp(), x$retrieved_at)
  x$plant_part = .plant_clean_context_value(x$plant_part)
  x$tissue = .plant_clean_context_value(x$tissue)
  x$method = .plant_clean_context_value(x$method)
  x$evidence_tier = .plant_normalize_evidence_tier(x$evidence_tier,
                                                   x$source_database,
                                                   x$matched_rank)
  x$confidence = .plant_normalize_confidence(x$confidence, x$evidence_tier)
  x$curation_flag = .uaf_squish_text(ifelse(is.na(x$curation_flag) |
                                              x$curation_flag == "",
                                            "unreviewed", x$curation_flag))
  x$plant_part_group = .plant_context_group(x$plant_part,
                                            x$plant_part_group,
                                            "plant_part")
  x$tissue_group = .plant_context_group(x$tissue, x$tissue_group,
                                        "tissue")
  x$method_group = .plant_method_group(x$method, x$method_group,
                                       x$source_database, x$evidence_tier)
  x$biological_context_status = .plant_context_status(x$plant_part_group,
                                                      x$tissue_group,
                                                      x$method_group)
  x$occurrence_status = .plant_occurrence_status(x$evidence_tier,
                                                 x$matched_rank,
                                                 x$source_database,
                                                 x$curation_flag)
  x$occurrence_basis = .plant_occurrence_basis(x$evidence_tier,
                                               x$matched_rank,
                                               x$source_database)
  x$evidence_quality_score = .plant_evidence_quality_score(
    x$occurrence_status, x$occurrence_basis, x$confidence,
    x$biological_context_status
  )
  x$analysis_ready = .plant_analysis_ready(x$occurrence_status,
                                           x$confidence,
                                           x$evidence_quality_score)
  x = x[!is.na(x$species) & x$species != "" &
          !is.na(x$compound_name) & x$compound_name != "", ,
        drop = FALSE]
  x = unique(x)
  row.names(x) = NULL
  x
}

.plant_normalize_literature = function(x, source_hint = "pubmed") {
  cols = .plant_literature_cols()
  if (!is.data.frame(x) || nrow(x) < 1) return(.uaf_empty_table(cols))
  names(x) = .plant_normalize_column_names(names(x))
  for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
  x = x[, cols, drop = FALSE]
  x$query_plant = .uaf_squish_text(x$query_plant)
  x$query_plant_clean = .plant_clean_name(x$query_plant)
  x$species = .uaf_squish_text(x$species)
  x$genus = ifelse(is.na(x$genus) | x$genus == "", .plant_genus(x$species),
                   .uaf_squish_text(x$genus))
  x$source_database = .uaf_squish_text(ifelse(is.na(x$source_database) |
                                                x$source_database == "",
                                              source_hint, x$source_database))
  x$retrieved_at = ifelse(is.na(x$retrieved_at) | x$retrieved_at == "",
                          .plant_timestamp(), x$retrieved_at)
  x$evidence_tier = .plant_normalize_evidence_tier(x$evidence_tier,
                                                   x$source_database,
                                                   "species")
  x$confidence = .plant_normalize_confidence(x$confidence, x$evidence_tier)
  x$curation_flag = .uaf_squish_text(ifelse(is.na(x$curation_flag) |
                                              x$curation_flag == "",
                                            "candidate", x$curation_flag))
  row.names(x) = NULL
  x
}

.plant_match_occurrences_to_queries = function(occurrences, plant_queries) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1) return(occurrences)
  lookup = plant_queries[, c("query_plant", "query_plant_clean", "species",
                             "genus", "family"), drop = FALSE]
  for (i in seq_len(nrow(occurrences))) {
    key = .plant_clean_name(occurrences$species[[i]])
    hit = lookup[lookup$query_plant_clean == key |
                   .plant_clean_name(lookup$species) == key, , drop = FALSE]
    if (nrow(hit) > 0) {
      occurrences$query_plant[[i]] = hit$query_plant[[1]]
      occurrences$query_plant_clean[[i]] = hit$query_plant_clean[[1]]
      occurrences$family[[i]] =
        .uaf_first_non_empty_text(occurrences$family[[i]], hit$family[[1]])
      occurrences$genus[[i]] =
        .uaf_first_non_empty_text(occurrences$genus[[i]], hit$genus[[1]])
    }
  }
  occurrences
}

.plant_match_literature_to_queries = function(literature, plant_queries) {
  literature = .plant_normalize_literature(literature)
  if (nrow(literature) < 1) return(literature)
  for (i in seq_len(nrow(literature))) {
    key = .plant_clean_name(literature$species[[i]])
    hit = plant_queries[plant_queries$query_plant_clean == key |
                          .plant_clean_name(plant_queries$species) == key,
                        , drop = FALSE]
    if (nrow(hit) > 0) {
      literature$query_plant[[i]] = hit$query_plant[[1]]
      literature$query_plant_clean[[i]] = hit$query_plant_clean[[1]]
      literature$family[[i]] =
        .uaf_first_non_empty_text(literature$family[[i]], hit$family[[1]])
      literature$genus[[i]] =
        .uaf_first_non_empty_text(literature$genus[[i]], hit$genus[[1]])
    }
  }
  literature
}

.plant_bind_occurrences = function(parts) {
  .plant_collapse_occurrence_evidence(
    .plant_normalize_occurrences(.plant_bind_tables(parts, .plant_occurrence_cols()))
  )
}

.plant_collapse_occurrence_evidence = function(occurrences) {
  if (!is.data.frame(occurrences) || nrow(occurrences) < 2) {
    return(occurrences)
  }
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  if (!all(dup_cols %in% names(occurrences))) return(occurrences)
  key_data = occurrences[dup_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(key_data, sep = "||"))
  groups = split(seq_len(nrow(occurrences)), key)
  rows = lapply(groups, function(idx) {
    if (length(idx) == 1) return(occurrences[idx, , drop = FALSE])
    .plant_merge_occurrence_group(occurrences[idx, , drop = FALSE])
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  .plant_normalize_occurrences(out)
}

.plant_merge_occurrence_group = function(group) {
  scores = suppressWarnings(as.numeric(group$evidence_quality_score))
  scores[!is.finite(scores)] = 0
  rank_score = .plant_occurrence_rank_score(group$occurrence_status) +
    .plant_confidence_score(group$confidence) + scores
  best = which.max(rank_score)
  row = group[best, , drop = FALSE]
  for (col in names(row)) {
    if (col %in% c("evidence_text", "curation_flag")) next
    if (length(.uaf_non_empty(row[[col]])) > 0) next
    values = .uaf_non_empty(group[[col]])
    if (length(values) > 0) {
      row[[col]] = values[[1]]
      if (identical(col, "plant_part") && "plant_part_group" %in% names(row)) {
        row$plant_part_group = NA_character_
      }
      if (identical(col, "tissue") && "tissue_group" %in% names(row)) {
        row$tissue_group = NA_character_
      }
      if (identical(col, "method") && "method_group" %in% names(row)) {
        row$method_group = NA_character_
      }
    }
  }
  row$evidence_text = .plant_truncate(
    .pubchem_collapse(unique(.uaf_non_empty(group$evidence_text))),
    1200
  )
  row$curation_flag = .pubchem_collapse(unique(.uaf_non_empty(group$curation_flag)))
  row
}

.plant_occurrence_rank_score = function(occurrence_status) {
  status = .plant_matrix_key(occurrence_status)
  out = rep(0, length(status))
  out[status == "curated_reported"] = 5
  out[status == "direct_reported"] = 4
  out[status == "literature_reported"] = 3
  out[status == "taxon_fallback"] = 2
  out[status == "candidate"] = 1
  out
}

.plant_bind_tables = function(parts, cols) {
  parts = parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) < 1) return(.uaf_empty_table(cols))
  normalized = lapply(parts, function(x) {
    if (!is.data.frame(x)) return(.uaf_empty_table(cols))
    for (col in cols) if (!col %in% names(x)) x[[col]] = NA_character_
    x[, cols, drop = FALSE]
  })
  out = do.call(rbind, normalized)
  if (is.null(out) || nrow(out) < 1) return(.uaf_empty_table(cols))
  row.names(out) = NULL
  out
}

.plant_compound_resolution = function(occurrences, categorate_result) {
  occurrences = .plant_normalize_occurrences(occurrences)
  compounds = sort(unique(.uaf_non_empty(occurrences$compound_name)))
  if (length(compounds) < 1) {
    return(.uaf_empty_table(.plant_compound_resolution_cols()))
  }
  props = .plant_categorate_properties(categorate_result)
  rows = lapply(compounds, function(compound) {
    clean = .plant_clean_compound(compound)
    hit = props[.plant_clean_compound(props$Query) == clean, , drop = FALSE]
    data.frame(
      compound_name = compound,
      compound_name_clean = clean,
      query_count = sum(occurrences$compound_name_clean == clean),
      resolved = nrow(hit) > 0 &&
        length(.uaf_non_empty(c(hit$CID, hit$InChIKey, hit$CanonicalSMILES,
                                hit$IsomericSMILES))) > 0,
      CID = .uaf_first_non_empty_text(hit$CID),
      InChIKey = .uaf_first_non_empty_text(hit$InChIKey),
      SMILES = .uaf_first_non_empty_text(hit$IsomericSMILES,
                                         hit$CanonicalSMILES),
      MolecularFormula = .uaf_first_non_empty_text(hit$MolecularFormula),
      resolution_source = ifelse(nrow(hit) > 0, "uafR_compound_enrichment",
                                 "unresolved"),
      notes = ifelse(nrow(hit) > 0, "", "No compound enrichment hit."),
      stringsAsFactors = FALSE
    )
  })
  out = do.call(rbind, rows)
  row.names(out) = NULL
  out
}

.plant_pubchem_only_enrichment = function(compounds, detail, cache, cache_dir,
                                          throttle, pubchem_fun = NULL,
                                          request_fun = NULL,
                                          batch_size = Inf,
                                          resume = TRUE,
                                          progress = FALSE,
                                          ...) {
  dots = list(...)
  batch_size = .plant_enrichment_batch_size(batch_size, length(compounds))
  if (batch_size < length(compounds)) {
    return(.plant_pubchem_only_enrichment_batched(
      compounds = compounds,
      detail = detail,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun,
      batch_size = batch_size,
      resume = resume,
      progress = progress,
      ...
    ))
  }
  pubchem_fun = pubchem_fun %||% pubchemProfile
  pubchem_profile = if (detail == "full") "full" else "safety"
  pubchem_sections = if (detail == "full") {
    c("Literature", "Patents")
  } else {
    c("Mass Spectrometry", "GC-MS", "MS-MS",
      "Names and Identifiers", "Chemical and Physical Properties",
      "Pharmacology and Biochemistry", "Drug and Medication Information")
  }
  pubchem_sources = c("LOTUS - the natural products occurrence database",
                      "Flavor and Extract Manufacturers Association (FEMA)",
                      "FDA/SPL Indexing Data",
                      "Medical Subject Headings (MeSH)")
  assay_detail_limit = dots$assay_detail_limit %||% 50
  trait_matrix_profile = dots$trait_matrix_profile %||% "core"
  trait_matrix_mode = dots$trait_matrix_mode %||% "binary"
  trait_matrix_min_confidence = dots$trait_matrix_min_confidence %||% 0
  trait_matrix_max_traits = dots$trait_matrix_max_traits %||% Inf

  pubchem_cache_dir = if (is.null(cache_dir)) NULL else
    file.path(cache_dir, "pubchem")
  pubchem = tryCatch(
    pubchem_fun(compounds = compounds,
                profile = pubchem_profile,
                sections = pubchem_sections,
                sources = pubchem_sources,
                cache = cache,
                cache_dir = pubchem_cache_dir,
                throttle = throttle,
                assay_detail_limit = assay_detail_limit,
                request_fun = request_fun),
    error = function(error) {
      warning("PubChem-only plant compound enrichment failed: ",
              conditionMessage(error), call. = FALSE)
      .categorate_empty_pubchem_profile(compounds, pubchem_profile)
    }
  )
  pubchem_profiles = .pubchem_extract_profiles(pubchem)
  kegg = .categorate_empty_kegg_profile()
  data_list = .plant_empty_categorate_data_list()
  normalized = .categorate_normalized_outputs(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    pubchem_profiles = pubchem_profiles,
    trait_matrix_profile = trait_matrix_profile,
    trait_matrix_mode = trait_matrix_mode,
    trait_matrix_min_confidence = trait_matrix_min_confidence,
    trait_matrix_max_traits = trait_matrix_max_traits
  )
  source_coverage = .categorate_source_coverage(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    pubchem_profiles = pubchem_profiles
  )
  derived_groups = .categorate_derived_groups(
    compounds = compounds,
    data_list = data_list,
    pubchem = pubchem,
    kegg = kegg,
    source_coverage = source_coverage,
    pubchem_profiles = pubchem_profiles,
    normalized = normalized
  )
  provenance = .categorate_research_provenance(pubchem = pubchem, kegg = kegg)
  result = list(
    PubChemIdentity = pubchem$identity,
    PubChemProperties = pubchem$properties,
    PubChemSynonyms = pubchem$synonyms,
    PubChemAnnotations = pubchem$annotations,
    PubChemSourceAnnotations = pubchem$source_annotations,
    PubChemSpectra = pubchem$spectra,
    PubChemSafety = pubchem$safety,
    PubChemExperimental = pubchem$experimental,
    PubChemBioactivity = pubchem$bioactivity,
    PubChemBioAssayDetails = pubchem$bioassay_details,
    PubChemIdentifiers = pubchem_profiles$PubChemIdentifierProfile,
    SafetyProfile = pubchem_profiles$SafetyProfile,
    FEMAProfile = pubchem_profiles$FEMAProfile,
    FDA_SPL_Profile = pubchem_profiles$FDA_SPL_Profile,
    LOTUSProfile = pubchem_profiles$LOTUSProfile,
    PubChemClassifications = pubchem_profiles$PubChemClassificationProfile,
    MeSHProfile = pubchem_profiles$MeSHProfile,
    LiteratureProfile = pubchem_profiles$LiteratureProfile,
    ChemicalTerms = normalized$ChemicalTerms,
    ChemicalTraits = normalized$ChemicalTraits,
    ChemicalTraitOntology = normalized$ChemicalTraitOntology,
    ChemicalTraitMatrix = normalized$ChemicalTraitMatrix,
    ChemicalTraitOntologyMatrix = normalized$ChemicalTraitOntologyMatrix,
    ChemicalTraitEvidence = normalized$ChemicalTraitEvidence,
    ChemicalTraitReport = normalized$ChemicalTraitReport,
    ChemicalTraitSummary = normalized$ChemicalTraitSummary,
    ChemicalTraitSimilarity = normalized$ChemicalTraitSimilarity,
    ChemicalClasses = normalized$ChemicalClasses,
    ChemicalMeasurements = normalized$ChemicalMeasurements,
    ChemicalMeasurementSummary = normalized$ChemicalMeasurementSummary,
    ChemicalHazards = normalized$ChemicalHazards,
    ChemicalUses = normalized$ChemicalUses,
    ChemicalBioassays = normalized$ChemicalBioassays,
    ChemicalBioactivities = normalized$ChemicalBioactivities,
    ChemicalTargets = normalized$ChemicalTargets,
    ChemicalPotencies = normalized$ChemicalPotencies,
    ChemicalTaxonomy = normalized$ChemicalTaxonomy,
    ChemicalOccurrences = normalized$ChemicalOccurrences,
    ChemicalPathwayRoles = normalized$ChemicalPathwayRoles,
    KEGGReactionParticipants = normalized$KEGGReactionParticipants,
    KEGGMatches = kegg$matches,
    KEGGRecords = kegg$records,
    KEGGIdentifiers = kegg$identifiers,
    KEGGPathways = kegg$pathways,
    KEGGReactions = kegg$reactions,
    KEGGEnzymes = kegg$enzymes,
    KEGGModules = kegg$modules,
    KEGGLinks = kegg$links,
    KEGGLinkMetadata = kegg$link_metadata,
    KEGGClassifications = kegg$classifications,
    SourceCoverage = source_coverage,
    DerivedGroups = derived_groups,
    Provenance = provenance,
    EnrichmentMode = "pubchem_only"
  )
  validation = validateCategorateResult(result)
  c(result, list(
    DataDictionary = validation$DataDictionary,
    TableQuality = validation$TableQuality,
    SourceDiagnostics = validation$SourceDiagnostics,
    ValidationIssues = validation$Issues,
    ValidationSummary = validation$Summary
  ))
}

.plant_enrichment_batch_size = function(batch_size, n_compounds) {
  batch_size = suppressWarnings(as.numeric(batch_size[[1]]))
  if (!is.finite(batch_size) || batch_size < 1) return(n_compounds)
  max(1, as.integer(batch_size))
}

.plant_pubchem_only_enrichment_batched = function(compounds, detail, cache,
                                                  cache_dir, throttle,
                                                  pubchem_fun, request_fun,
                                                  batch_size, resume,
                                                  progress, ...) {
  groups = split(compounds, ceiling(seq_along(compounds) / batch_size))
  results = vector("list", length(groups))
  batch_dir = if (is.null(cache_dir)) {
    file.path(.plant_cache_dir(NULL), "compound_enrichment",
              "pubchem_only_batches")
  } else {
    file.path(cache_dir, "pubchem_only_batches")
  }
  if (isTRUE(cache)) dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  for (i in seq_along(groups)) {
    batch = groups[[i]]
    cache_file = file.path(
      batch_dir,
      paste0("batch_", sprintf("%04d", i), "_",
             .pubchem_url_hash(paste(c(detail, batch), collapse = "\r")),
             ".rds")
    )
    if (isTRUE(cache) && isTRUE(resume) && file.exists(cache_file)) {
      if (isTRUE(progress)) {
        message("uafR plant enrichment batch ", i, "/", length(groups),
                ": using cached result")
      }
      results[[i]] = readRDS(cache_file)
      next
    }
    if (isTRUE(progress)) {
      message("uafR plant enrichment batch ", i, "/", length(groups),
              ": resolving ", length(batch), " compound(s)")
    }
    result = .plant_pubchem_only_enrichment(
      compounds = batch,
      detail = detail,
      cache = cache,
      cache_dir = cache_dir,
      throttle = throttle,
      pubchem_fun = pubchem_fun,
      request_fun = request_fun,
      batch_size = Inf,
      resume = FALSE,
      progress = FALSE,
      ...
    )
    if (isTRUE(cache)) saveRDS(result, cache_file)
    results[[i]] = result
  }
  .plant_merge_categorate_results(results, compounds,
                                  enrichment_mode = "pubchem_only_batched")
}

.plant_merge_categorate_results = function(results, compounds,
                                           enrichment_mode) {
  results = results[vapply(results, is.list, logical(1))]
  if (length(results) < 1) return(NULL)
  names_all = unique(unlist(lapply(results, names), use.names = FALSE))
  merged = list()
  for (name in names_all) {
    values = lapply(results, `[[`, name)
    values = values[!vapply(values, is.null, logical(1))]
    if (length(values) < 1) next
    if (all(vapply(values, is.data.frame, logical(1)))) {
      merged[[name]] = .plant_bind_flexible_tables(values)
    } else {
      merged[[name]] = values[[1]]
    }
  }
  merged$EnrichmentMode = enrichment_mode
  merged$BatchCount = length(results)
  merged$BatchCompoundCount = length(compounds)
  validation = validateCategorateResult(merged)
  c(merged, list(
    DataDictionary = validation$DataDictionary,
    TableQuality = validation$TableQuality,
    SourceDiagnostics = validation$SourceDiagnostics,
    ValidationIssues = validation$Issues,
    ValidationSummary = validation$Summary
  ))
}

.plant_bind_flexible_tables = function(parts) {
  parts = parts[vapply(parts, is.data.frame, logical(1))]
  if (length(parts) < 1) return(data.frame())
  cols = unique(unlist(lapply(parts, names), use.names = FALSE))
  normalized = lapply(parts, function(part) {
    for (col in cols) if (!col %in% names(part)) part[[col]] = NA_character_
    part[, cols, drop = FALSE]
  })
  out = do.call(rbind, normalized)
  out = unique(out)
  row.names(out) = NULL
  out
}

.plant_empty_categorate_data_list = function() {
  list(
    reactives = .uaf_empty_table(c("Chemical", "reactives")),
    LOTUS = .uaf_empty_table(c("Chemical", "LOTUS")),
    KEGG = .uaf_empty_table(c("Chemical", "KEGG")),
    FEMA = .uaf_empty_table(c("Chemical", "FEMA")),
    FDA_SPL = .uaf_empty_table(c("Chemical", "FDA_SPL"))
  )
}

.plant_enrichment_source = function(categorate_result) {
  mode = if (is.list(categorate_result)) {
    .uaf_first_non_empty_text(categorate_result$EnrichmentMode)
  } else {
    NA_character_
  }
  if (!is.na(mode) && mode %in% c("pubchem_only", "pubchem_only_batched")) {
    return("pubchemProfile")
  }
  "categorate"
}

.plant_enrichment_note = function(categorate_result) {
  mode = if (is.list(categorate_result)) {
    .uaf_first_non_empty_text(categorate_result$EnrichmentMode)
  } else {
    NA_character_
  }
  if (identical(mode, "pubchem_only_batched")) {
    return("PubChem-only compound enrichment completed in resumable batches.")
  }
  if (identical(mode, "pubchem_only")) {
    return("PubChem-only compound enrichment completed.")
  }
  "Compound enrichment completed or supplied."
}

.plant_categorate_properties = function(categorate_result) {
  cols = c("Query", "CID", "MolecularFormula", "InChIKey",
           "CanonicalSMILES", "IsomericSMILES")
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$PubChemProperties)) {
    props = categorate_result$PubChemProperties
  } else if (is.list(categorate_result) &&
             is.data.frame(categorate_result$properties)) {
    props = categorate_result$properties
  } else {
    return(.uaf_empty_table(cols))
  }
  for (col in cols) if (!col %in% names(props)) props[[col]] = NA_character_
  props[, cols, drop = FALSE]
}

.plant_trait_evidence = function(occurrences, categorate_result) {
  occurrences = .plant_normalize_occurrences(occurrences)
  if (nrow(occurrences) < 1 || !is.list(categorate_result) ||
      !is.data.frame(categorate_result$ChemicalTraitEvidence)) {
    return(.uaf_empty_table(.plant_trait_evidence_cols()))
  }
  evidence = categorate_result$ChemicalTraitEvidence
  if (!"Query" %in% names(evidence)) {
    return(.uaf_empty_table(.plant_trait_evidence_cols()))
  }
  rows = list()
  for (i in seq_len(nrow(occurrences))) {
    occ = occurrences[i, , drop = FALSE]
    hit = evidence[.plant_clean_compound(evidence$Query) ==
                     occ$compound_name_clean, , drop = FALSE]
    if (nrow(hit) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      species = occ$species,
      species_slug = .plant_slug(occ$species),
      compound_name = occ$compound_name,
      compound_name_clean = occ$compound_name_clean,
      EvidenceType = .plant_col_or_default(hit, "EvidenceType", NA_character_),
      AnalysisKey = .plant_col_or_default(hit, "AnalysisKey", NA_character_),
      SourceDatabase = .plant_col_or_default(hit, "SourceDatabase", NA_character_),
      EvidenceText = .plant_col_or_default(hit, "EvidenceText", NA_character_),
      EvidenceURL = .plant_col_or_default(hit, "EvidenceURL", NA_character_),
      Confidence = .plant_col_or_default(hit, "Confidence", NA_character_),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, .plant_trait_evidence_cols())
}

.plant_summary_trait_info = function(group, categorate_result) {
  compounds = unique(.uaf_non_empty(group$compound_name_clean))
  traits = if (is.list(categorate_result) &&
               is.data.frame(categorate_result$ChemicalTraits)) {
    categorate_result$ChemicalTraits
  } else {
    .uaf_empty_table(c("Query", "TraitType", "TraitGroup", "TraitValue",
                       "Confidence", "SourceDatabase"))
  }
  classes = if (is.list(categorate_result) &&
                is.data.frame(categorate_result$ChemicalClasses)) {
    categorate_result$ChemicalClasses
  } else {
    .uaf_empty_table(c("Query", "ClassName", "ClassGroup",
                       "Superclass", "Class", "Subclass"))
  }
  keep_traits = if ("Query" %in% names(traits)) {
    .plant_clean_compound(traits$Query) %in% compounds
  } else logical(nrow(traits))
  keep_classes = if ("Query" %in% names(classes)) {
    .plant_clean_compound(classes$Query) %in% compounds
  } else logical(nrow(classes))
  traits = traits[keep_traits, , drop = FALSE]
  classes = classes[keep_classes, , drop = FALSE]
  text = tolower(.pubchem_collapse(c(group$compound_name, group$evidence_text,
                                     traits$TraitType, traits$TraitGroup,
                                     traits$TraitValue, classes$ClassName,
                                     classes$ClassGroup)))
  fraction = function(pattern) {
    if (length(compounds) < 1) return(0)
    hits = vapply(compounds, function(compound) {
      compound_text = tolower(.pubchem_collapse(c(
        group$compound_name[group$compound_name_clean == compound],
        traits$TraitValue[.plant_clean_compound(traits$Query) == compound],
        traits$TraitGroup[.plant_clean_compound(traits$Query) == compound]
      )))
      grepl(pattern, compound_text, perl = TRUE)
    }, logical(1))
    mean(hits)
  }
  list(
    superclasses = .pubchem_collapse(.plant_col_or_default(classes, "Superclass", NA_character_)),
    classes = .pubchem_collapse(c(.plant_col_or_default(classes, "Class", NA_character_),
                                  .plant_col_or_default(classes, "ClassName", NA_character_))),
    subclasses = .pubchem_collapse(.plant_col_or_default(classes, "Subclass", NA_character_)),
    dominant_classes = .plant_top_terms(c(.plant_col_or_default(classes, "ClassName", NA_character_),
                                          .plant_col_or_default(traits, "TraitGroup", NA_character_))),
    kingdoms = ifelse(grepl("plantae|plant", text), "Plantae", NA_character_),
    volatile_fraction = fraction("volatile|odor|aroma|essential oil|gc-ms|gc ms"),
    lipophilic_fraction = fraction("lipophilic|xlogp|logp|terpene|fatty|alkane|lipid"),
    oxygenated_fraction = fraction("oxygen|hydroxy|phenol|acid|ester|alcohol|ketone|aldehyde"),
    nitrogenous_fraction = fraction("nitrogen|alkaloid|amine|amide|nitrile"),
    sulfur_fraction = fraction("sulfur|thiol|sulfide|sulfoxide"),
    halogen_fraction = fraction("chloro|bromo|fluoro|iodo|halogen"),
    kegg_pathway_groups = .pubchem_collapse(.plant_col_or_default(traits, "PathwayGroup", NA_character_)),
    chemical_trait_count = nrow(traits),
    high_confidence_trait_count = sum(tolower(.plant_col_or_default(traits, "Confidence", "")) == "high")
  )
}

.plant_wide_matrix = function(long, id_col, mode, max_traits) {
  if (nrow(long) < 1) return(.uaf_empty_table(c(id_col)))
  traits = names(sort(table(long$trait), decreasing = TRUE))
  max_traits = suppressWarnings(as.numeric(max_traits[[1]]))
  if (is.finite(max_traits)) traits = utils::head(traits, max_traits)
  ids = sort(unique(.uaf_non_empty(long$id)))
  out = data.frame(stats::setNames(list(ids), id_col), stringsAsFactors = FALSE)
  for (trait in traits) {
    values = vapply(ids, function(id) {
      vals = long$value[long$id == id & long$trait == trait]
      if (mode == "binary") as.numeric(length(vals) > 0)
      else if (mode == "count") length(vals)
      else max(vals, na.rm = TRUE)
    }, numeric(1))
    values[!is.finite(values)] = 0
    out[[trait]] = values
  }
  row.names(out) = NULL
  out
}

.plant_chemical_trait_matrix_long = function(occurrences, categorate_result,
                                             level, profile, min_confidence) {
  cols = c("id", "trait", "value")
  if (!is.list(categorate_result) ||
      !is.data.frame(categorate_result$ChemicalTraits) ||
      nrow(categorate_result$ChemicalTraits) < 1 ||
      nrow(occurrences) < 1) {
    return(.uaf_empty_table(cols))
  }
  traits = .normalized_prepare_traits(categorate_result$ChemicalTraits)
  traits = .normalized_filter_trait_matrix_source(
    chemical_traits = traits,
    profile = .plant_chemical_trait_profile(profile),
    min_confidence = .normalized_trait_confidence_threshold(min_confidence)
  )
  if (nrow(traits) < 1) return(.uaf_empty_table(cols))
  traits$query_clean = .plant_clean_compound(traits$Query)
  rows = list()
  for (i in seq_len(nrow(occurrences))) {
    occurrence = occurrences[i, , drop = FALSE]
    id = .uaf_first_non_empty_text(occurrence[[level]],
                                   occurrence$species)
    if (is.na(id) || id == "") next
    hits = traits[traits$query_clean == occurrence$compound_name_clean, ,
                  drop = FALSE]
    if (nrow(hits) < 1) next
    rows[[length(rows) + 1]] = data.frame(
      id = id,
      trait = paste0("chemtrait__", hits$MatrixKey),
      compound_name_clean = occurrence$compound_name_clean,
      value = suppressWarnings(as.numeric(hits$ConfidenceScore)),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(cols))
  out = do.call(rbind, rows)
  out$value[is.na(out$value)] = 0
  out = unique(out[, c("id", "trait", "compound_name_clean", "value"),
                   drop = FALSE])
  out = out[, cols, drop = FALSE]
  row.names(out) = NULL
  out
}

.plant_context_matrix_long = function(occurrences, level, group_col, prefix,
                                      exclude = c("unknown")) {
  cols = c("id", "trait", "value")
  if (!is.data.frame(occurrences) || nrow(occurrences) < 1 ||
      !group_col %in% names(occurrences)) {
    return(.uaf_empty_table(cols))
  }
  keep = !is.na(occurrences[[group_col]]) &
    occurrences[[group_col]] != "" &
    !occurrences[[group_col]] %in% exclude
  rows = occurrences[keep, , drop = FALSE]
  if (nrow(rows) < 1) return(.uaf_empty_table(cols))
  id = rows[[level]]
  id[is.na(id) | id == ""] = rows$species[is.na(id) | id == ""]
  value = suppressWarnings(as.numeric(rows$evidence_quality_score))
  missing_value = !is.finite(value)
  value[missing_value] = .plant_confidence_score(rows$confidence[missing_value])
  out = data.frame(
    id = id,
    trait = paste0(prefix, "__", .plant_matrix_key(rows[[group_col]])),
    value = value,
    stringsAsFactors = FALSE
  )
  out = unique(out)
  row.names(out) = NULL
  out
}

.plant_chemical_trait_profile = function(profile) {
  if (identical(profile, "metabolism")) return("kegg")
  profile
}

.plant_export_tables = function(x, tables, include_empty, max_cell_chars) {
  default = c("PlantQueries", "PlantNameResolution", "ProviderDiagnostics",
              "PlantCompoundOccurrences", "LiteratureCandidates",
              "CompoundResolution", "SpeciesChemistrySummary",
              "SpeciesChemistryMatrix", "TraitEvidence", "ValidationSummary",
              "ValidationIssues", "DataDictionary", "Provenance")
  available = list(
    PlantQueries = x$PlantQueries,
    PlantNameResolution = x$PlantNameResolution,
    ProviderDiagnostics = x$ProviderDiagnostics,
    PlantCompoundOccurrences = x$PlantCompoundOccurrences,
    LiteratureCandidates = x$LiteratureCandidates,
    CompoundResolution = x$CompoundResolution,
    SpeciesChemistrySummary = x$SpeciesChemistrySummary,
    SpeciesChemistryMatrix = x$SpeciesChemistryMatrix,
    TraitEvidence = x$TraitEvidence,
    ValidationSummary = if (is.list(x$Validation)) x$Validation$Summary else NULL,
    ValidationIssues = if (is.list(x$Validation)) x$Validation$Issues else NULL,
    DataDictionary = x$DataDictionary,
    Provenance = x$Provenance
  )
  requested = .uaf_non_empty(tables)
  table_names = if (length(requested) > 0) requested else default
  out = list()
  for (name in table_names) {
    table = available[[name]]
    if (!is.data.frame(table)) next
    if (!isTRUE(include_empty) && nrow(table) < 1) next
    out[[name]] = .categorate_export_prepare_table(table, max_cell_chars)
  }
  out
}

.plant_apply_optional_filter = function(keep, values, allowed) {
  allowed = .uaf_non_empty(allowed)
  if (length(allowed) < 1 || "all" %in% tolower(allowed)) return(keep)
  keep & values %in% allowed
}

.plant_review_cols = function() {
  c("review_id", "review_decision", "reviewed_by", "reviewed_at",
    "review_note", "proposed_evidence_tier", "proposed_confidence",
    "proposed_source_database", "proposed_citation_or_url",
    "proposed_plant_part", "proposed_tissue", "proposed_method",
    "review_key", .plant_occurrence_cols())
}

.plant_occurrence_key = function(occurrences) {
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  key_data = occurrences[dup_cols]
  key_data[] = lapply(key_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  do.call(paste, c(key_data, sep = "||"))
}

.plant_review_match_indices = function(occurrences, review) {
  keys = .plant_occurrence_key(occurrences)
  review_key = if ("review_key" %in% names(review)) {
    as.character(review$review_key[[1]])
  } else {
    NA_character_
  }
  if (is.na(review_key) || review_key == "") review_key = NA_character_
  if (!is.na(review_key)) return(which(keys == review_key))
  dup_cols = .plant_duplicate_key_cols("PlantCompoundOccurrences")
  if (!all(dup_cols %in% names(review))) return(integer())
  review_data = review[dup_cols]
  review_data[] = lapply(review_data, function(value) {
    value = .uaf_squish_text(value)
    value[is.na(value)] = ""
    value
  })
  key = do.call(paste, c(review_data, sep = "||"))
  which(keys == key)
}

.plant_review_decision = function(x) {
  decision = tolower(.uaf_first_non_empty_text(x))
  if (is.na(decision)) return("needs_review")
  decision = gsub("[^a-z0-9]+", "_", decision)
  decision = gsub("^_|_$", "", decision)
  if (decision %in% c("accept", "promote", "promote_reported",
                      "accept_curated", "curated")) {
    return("promote_curated")
  }
  if (decision %in% c("remove", "drop")) return("reject")
  if (decision %in% c("update", "context", "edit_context")) {
    return("update_context")
  }
  if (decision %in% c("candidate", "keep_as_candidate")) {
    return("keep_candidate")
  }
  decision
}

.plant_apply_review_context = function(occurrences, idx, review) {
  update = function(col, value) {
    value = .uaf_first_non_empty_text(value)
    if (!is.na(value) && col %in% names(occurrences)) {
      occurrences[[col]][idx] <<- value
      TRUE
    } else {
      FALSE
    }
  }
  if (update("plant_part", review$proposed_plant_part) &&
      "plant_part_group" %in% names(occurrences)) {
    occurrences$plant_part_group[idx] = NA_character_
  }
  if (update("tissue", review$proposed_tissue) &&
      "tissue_group" %in% names(occurrences)) {
    occurrences$tissue_group[idx] = NA_character_
  }
  if (update("method", review$proposed_method) &&
      "method_group" %in% names(occurrences)) {
    occurrences$method_group[idx] = NA_character_
  }
  update("confidence", review$proposed_confidence)
  update("source_database", review$proposed_source_database)
  update("evidence_url", review$proposed_citation_or_url)
  update("evidence_tier", review$proposed_evidence_tier)
  occurrences
}

.plant_rebuild_filtered_result = function(x, occurrences) {
  occurrences = .plant_normalize_occurrences(occurrences)
  out = x
  out$PlantCompoundOccurrences = occurrences
  out$CompoundResolution = .plant_compound_resolution(occurrences,
                                                      x$CategorateResult)
  out$TraitEvidence = .plant_trait_evidence(occurrences, x$CategorateResult)
  out$SpeciesChemistrySummary = summarizePlantPhytochemistry(
    plant_compounds = occurrences,
    categorate_result = x$CategorateResult,
    compound_resolution = out$CompoundResolution,
    plant_queries = x$PlantQueries,
    provider_diagnostics = x$ProviderDiagnostics
  )
  out$SpeciesChemistryMatrix = plantPhytochemistryMatrix(
    list(PlantCompoundOccurrences = occurrences,
         CategorateResult = x$CategorateResult),
    level = "species",
    profile = "core",
    mode = "binary",
    min_confidence = "medium"
  )
  out$Validation = validatePlantPhytochemistryResult(out)
  class(out) = unique(c("uaf_plant_phytochemistry", class(x)))
  out
}

.plant_validation_result = function(table_quality, provider_diagnostics,
                                    issues, dictionary) {
  summary = .plant_validation_summary(table_quality, issues)
  out = list(
    Summary = summary,
    TableQuality = table_quality,
    ProviderDiagnostics = provider_diagnostics,
    Issues = issues,
    DataDictionary = dictionary
  )
  class(out) = c("uaf_plant_phytochemistry_validation", class(out))
  out
}

.plant_validation_summary = function(table_quality, issues) {
  errors = sum(issues$severity == "error")
  warnings = sum(issues$severity == "warning")
  data.frame(
    Status = ifelse(errors > 0, "fail", ifelse(warnings > 0, "warning", "pass")),
    TablesExpected = nrow(table_quality),
    TablesPresent = sum(table_quality$Present == "Yes"),
    IssueCount = nrow(issues),
    ErrorCount = errors,
    WarningCount = warnings,
    stringsAsFactors = FALSE
  )
}

.plant_table_quality = function(table_name, table, dictionary, strict) {
  if (table_name == "Validation" && is.list(table)) {
    missing_required = setdiff(dictionary$Column[dictionary$Required],
                               names(table))
    return(data.frame(
      Table = table_name,
      Present = "Yes",
      RowCount = NA_integer_,
      ColumnCount = length(table),
      MissingRequiredColumns = .pubchem_collapse(missing_required),
      MissingOptionalColumns = "",
      DuplicateKeyColumns = "",
      Status = ifelse(length(missing_required) > 0, "fail", "pass"),
      stringsAsFactors = FALSE
    ))
  }
  present = is.data.frame(table)
  required = dictionary$Column[dictionary$Required]
  missing_required = if (present) setdiff(required, names(table)) else required
  optional = dictionary$Column[!dictionary$Required]
  missing_optional = if (present && isTRUE(strict)) {
    setdiff(optional, names(table))
  } else character()
  data.frame(
    Table = table_name,
    Present = .uaf_yes_no(present),
    RowCount = if (present) nrow(table) else NA_integer_,
    ColumnCount = if (present) ncol(table) else NA_integer_,
    MissingRequiredColumns = .pubchem_collapse(missing_required),
    MissingOptionalColumns = .pubchem_collapse(missing_optional),
    DuplicateKeyColumns = .pubchem_collapse(.plant_duplicate_key_cols(table_name)),
    Status = ifelse(!present || length(missing_required) > 0, "fail", "pass"),
    stringsAsFactors = FALSE
  )
}

.plant_validate_table = function(table_name, table, dictionary, strict) {
  rows = list()
  if (table_name == "Validation" && is.list(table)) {
    missing_required = setdiff(dictionary$Column[dictionary$Required],
                               names(table))
    for (col in missing_required) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("error", table_name, col,
                                "Required validation component is missing",
                                "present", "absent", NA_integer_,
                                NA_character_)
    }
    return(.plant_bind_tables(rows, .plant_validation_issue_cols()))
  }
  if (!is.data.frame(table)) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("error", table_name, NA_character_,
                              "Expected table is missing", "data.frame",
                              class(table)[[1]], NA_integer_, NA_character_)
    return(.plant_bind_tables(rows, .plant_validation_issue_cols()))
  }
  required = dictionary$Column[dictionary$Required]
  missing_required = setdiff(required, names(table))
  for (col in missing_required) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("error", table_name, col,
                              "Required column is missing", "present",
                              "absent", nrow(table), NA_character_)
  }
  if (isTRUE(strict)) {
    missing_optional = setdiff(dictionary$Column[!dictionary$Required],
                               names(table))
    for (col in missing_optional) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name, col,
                                "Optional documented column is missing",
                                "present", "absent", nrow(table),
                                NA_character_)
    }
  }
  present_dictionary = dictionary[dictionary$Column %in% names(table), ,
                                  drop = FALSE]
  for (i in seq_len(nrow(present_dictionary))) {
    allowed = .uaf_non_empty(strsplit(
      .uaf_first_non_empty_text(present_dictionary$AllowedValues[[i]]),
      ";", fixed = TRUE
    )[[1]])
    if (length(allowed) < 1) next
    col = present_dictionary$Column[[i]]
    bad = setdiff(unique(.uaf_non_empty(table[[col]])), allowed)
    if (length(bad) > 0) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name, col,
                                "Column contains values outside documented set",
                                .pubchem_collapse(allowed),
                                .pubchem_collapse(utils::head(bad, 8)),
                                nrow(table), NA_character_)
    }
  }
  dup_cols = .plant_duplicate_key_cols(table_name)
  if (length(dup_cols) > 0 && all(dup_cols %in% names(table)) &&
      nrow(table) > 0) {
    key = do.call(paste, c(table[dup_cols], sep = "\r"))
    key = key[!is.na(key) & key != ""]
    dup = unique(key[duplicated(key)])
    if (length(dup) > 0) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", table_name,
                                .pubchem_collapse(dup_cols),
                                "Duplicate evidence keys were found",
                                "unique rows",
                                paste(length(dup), "duplicate key(s)"),
                                nrow(table),
                                .pubchem_collapse(utils::head(dup, 5)))
    }
  }
  .plant_bind_tables(rows, .plant_validation_issue_cols())
}

.plant_extra_validation_issues = function(x, issues) {
  rows = list(issues)
  if (is.data.frame(x$PlantNameResolution) &&
      any(x$PlantNameResolution$query_status == "name_needs_review")) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "PlantNameResolution",
                              "query_status",
                              "Some plant names need taxonomic review",
                              "parsed_species",
                              "name_needs_review",
                              nrow(x$PlantNameResolution),
                              .pubchem_collapse(x$PlantNameResolution$query_plant[
                                x$PlantNameResolution$query_status ==
                                  "name_needs_review"
                              ]))
  }
  if (is.data.frame(x$CompoundResolution) &&
      any(!(x$CompoundResolution$resolved %in% TRUE))) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "CompoundResolution", "resolved",
                              "Some compounds are unresolved",
                              "all TRUE",
                              paste(sum(!(x$CompoundResolution$resolved %in% TRUE)),
                                    "unresolved"),
                              nrow(x$CompoundResolution),
                              .pubchem_collapse(x$CompoundResolution$compound_name[
                                !(x$CompoundResolution$resolved %in% TRUE)
                              ]))
  }
  if (is.data.frame(x$PlantCompoundOccurrences) &&
      "evidence_quality_score" %in% names(x$PlantCompoundOccurrences)) {
    score = suppressWarnings(as.numeric(
      x$PlantCompoundOccurrences$evidence_quality_score
    ))
    bad_score = is.na(score) | score < 0 | score > 1
    if (any(bad_score)) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", "PlantCompoundOccurrences",
                                "evidence_quality_score",
                                "Evidence quality scores must be between 0 and 1",
                                "0 <= score <= 1",
                                .pubchem_collapse(utils::head(score[bad_score], 8)),
                                nrow(x$PlantCompoundOccurrences),
                                NA_character_)
    }
    bad_ready = x$PlantCompoundOccurrences$analysis_ready == "Yes" &
      !x$PlantCompoundOccurrences$occurrence_status %in%
      c("direct_reported", "curated_reported")
    if (any(bad_ready, na.rm = TRUE)) {
      rows[[length(rows) + 1]] =
        .plant_validation_issue("warning", "PlantCompoundOccurrences",
                                "analysis_ready",
                                "Analysis-ready rows should be direct or curated reported occurrences",
                                "direct_reported or curated_reported",
                                .pubchem_collapse(unique(
                                  x$PlantCompoundOccurrences$occurrence_status[bad_ready]
                                )),
                                nrow(x$PlantCompoundOccurrences),
                                .pubchem_collapse(utils::head(
                                  x$PlantCompoundOccurrences$compound_name[bad_ready],
                                  8
                                )))
    }
  }
  if (is.data.frame(x$ProviderDiagnostics) &&
      any(x$ProviderDiagnostics$status %in% c("not_implemented",
                                              "not_queried",
                                              "warning"))) {
    rows[[length(rows) + 1]] =
      .plant_validation_issue("warning", "ProviderDiagnostics", "status",
                              "One or more providers did not return usable records",
                              "ok or no_records",
                              .pubchem_collapse(unique(x$ProviderDiagnostics$status)),
                              nrow(x$ProviderDiagnostics),
                              .pubchem_collapse(x$ProviderDiagnostics$provider[
                                x$ProviderDiagnostics$status %in%
                                  c("not_implemented", "not_queried",
                                    "warning")
                              ]))
  }
  .plant_bind_tables(rows, .plant_validation_issue_cols())
}

.plant_validation_issue = function(severity, table, column, issue, expected,
                                   observed, row_count, examples) {
  data.frame(
    severity = severity,
    table = table,
    column = column,
    issue = issue,
    expected = expected,
    observed = observed,
    row_count = row_count,
    examples = examples,
    stringsAsFactors = FALSE
  )
}

.plant_duplicate_key_cols = function(table_name) {
  switch(table_name,
         PlantCompoundOccurrences = c("species", "compound_name_clean",
                                      "source_database", "source_record_id",
                                      "pmid", "evidence_url"),
         LiteratureCandidates = c("species", "source_database", "pmid",
                                  "chemical_mention"),
         character())
}

.plant_empty_occurrences = function() .uaf_empty_table(.plant_occurrence_cols())

.plant_provenance = function(step, source, query, source_url, record_count,
                             notes, cache_file = NA_character_) {
  data.frame(
    step = step,
    source = source,
    query = query,
    source_url = source_url,
    cache_file = cache_file,
    retrieved_at = .plant_timestamp(),
    record_count = suppressWarnings(as.integer(record_count)),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

.plant_pubmed_query = function(species) {
  paste0('"', species, '"[Title/Abstract] AND (phytochemical* OR metabolite* OR ',
         '"secondary metabolite*" OR "chemical profile" OR "metabolic profile" OR ',
         '"natural product*" OR GC-MS OR "GC MS" OR LC-MS OR "LC MS" OR HPLC OR ',
         '"essential oil" OR volatile* OR flavonoid* OR alkaloid* OR terpene* OR ',
         'phenolic* OR "root exudate*")')
}

.plant_ncbi_url = function(endpoint, params, api_key = "") {
  params = params[!is.na(params) & params != ""]
  if (length(.uaf_non_empty(api_key)) > 0) {
    params = c(params, api_key = api_key)
  }
  paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/", endpoint, "?",
         paste(names(params), utils::URLencode(params, reserved = TRUE),
               sep = "=", collapse = "&"))
}

.plant_fetch_json = function(url, cache, cache_dir, throttle, request_fun) {
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".json"))
  if (isTRUE(cache) && file.exists(cache_file)) {
    txt = paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n")
    return(jsonlite::fromJSON(txt, simplifyVector = FALSE))
  }
  result = if (is.null(request_fun)) {
    con = base::url(url, open = "rb")
    on.exit(close(con), add = TRUE)
    paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    request_fun(url)
  }
  if (is.character(result)) {
    txt = paste(result, collapse = "\n")
    if (isTRUE(cache)) {
      dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
      writeLines(txt, cache_file, useBytes = TRUE)
    }
    if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
  } else {
    result
  }
}

.plant_fetch_text = function(url, cache, cache_dir, throttle, request_fun) {
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".txt"))
  if (isTRUE(cache) && file.exists(cache_file)) {
    return(paste(readLines(cache_file, warn = FALSE, encoding = "UTF-8"),
                 collapse = "\n"))
  }
  result = if (is.null(request_fun)) {
    con = base::url(url, open = "rb")
    on.exit(close(con), add = TRUE)
    paste(readLines(con, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    request_fun(url)
  }
  txt = if (is.list(result)) {
    jsonlite::toJSON(result, auto_unbox = TRUE)
  } else {
    paste(as.character(result), collapse = "\n")
  }
  if (isTRUE(cache)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(txt, cache_file, useBytes = TRUE)
  }
  if (!is.na(throttle) && throttle > 0) Sys.sleep(throttle)
  txt
}

.plant_pubmed_ids = function(search) {
  ids = tryCatch(search$esearchresult$idlist, error = function(error) NULL)
  .uaf_non_empty(unlist(ids, use.names = FALSE))
}

.plant_pubmed_summary_rows = function(query_row, ids, summary) {
  rows = list()
  for (id in ids) {
    item = tryCatch(summary$result[[id]], error = function(error) NULL)
    title = .uaf_first_non_empty_text(item$title)
    doi = .plant_article_id(item, "doi")
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      source_database = "PubMed",
      source_record_id = id,
      pmid = id,
      doi = doi,
      title = title,
      abstract = NA_character_,
      chemical_mention = NA_character_,
      species_mention = query_row$species,
      evidence_text = title,
      evidence_url = paste0("https://pubmed.ncbi.nlm.nih.gov/", id, "/"),
      retrieved_at = .plant_timestamp(),
      confidence = "low",
      curation_flag = "literature_candidate",
      evidence_tier = "direct_species_literature",
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_tables(rows, .plant_literature_cols())
}

.plant_taxon_query_terms = function(query_row) {
  fallbacks = .uaf_non_empty(strsplit(
    .uaf_first_non_empty_text(query_row$taxon_fallback),
    ";",
    fixed = TRUE
  )[[1]])
  fallbacks = tolower(.uaf_squish_text(fallbacks))
  ranks = unique(c("species", fallbacks))
  rows = list()
  if ("species" %in% ranks) {
    rows[[length(rows) + 1]] = data.frame(
      term = .uaf_first_non_empty_text(query_row$species,
                                       query_row$query_plant),
      rank = "species",
      stringsAsFactors = FALSE
    )
  }
  if ("genus" %in% ranks && length(.uaf_non_empty(query_row$genus)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      term = .uaf_first_non_empty_text(query_row$genus),
      rank = "genus",
      stringsAsFactors = FALSE
    )
  }
  if ("family" %in% ranks && length(.uaf_non_empty(query_row$family)) > 0) {
    rows[[length(rows) + 1]] = data.frame(
      term = .uaf_first_non_empty_text(query_row$family),
      rank = "family",
      stringsAsFactors = FALSE
    )
  }
  out = .plant_bind_tables(rows, c("term", "rank"))
  out = out[!is.na(out$term) & out$term != "", , drop = FALSE]
  unique(out)
}

.plant_lotus_rows = function(query_row, result, url, max_records) {
  items = .plant_json_records(result)
  if (length(items) < 1) return(.plant_empty_occurrences())
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records)) items = utils::head(items, max_records)
  rows = list()
  for (item in items) {
    flat = .plant_flatten_record(item)
    if (length(flat) < 1) next
    evidence_text = .plant_truncate(.pubchem_collapse(flat), 800)
    taxon_values = flat[grepl("organism|taxon|species|genus|family",
                              names(flat), ignore.case = TRUE)]
    taxon_text = .pubchem_collapse(taxon_values)
    matched_rank = .plant_matched_rank_from_text(query_row, taxon_text)
    if (is.na(matched_rank)) {
      matched_rank = .plant_matched_rank_from_text(query_row, evidence_text)
    }
    if (is.na(matched_rank)) next
    compound_name = .plant_lotus_compound_name(flat, query_row)
    if (is.na(compound_name)) next
    record_id = .plant_lotus_record_id(flat)
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      matched_taxon = ifelse(matched_rank == "genus", query_row$genus,
                             ifelse(matched_rank == "family",
                                    query_row$family, query_row$species)),
      matched_rank = matched_rank,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      compound_name = compound_name,
      compound_name_clean = .plant_clean_compound(compound_name),
      compound_id = record_id,
      compound_id_type = ifelse(is.na(record_id), NA_character_,
                                "LOTUS_or_Wikidata"),
      source_database = "LOTUS",
      source_record_id = record_id,
      evidence_text = evidence_text,
      evidence_url = url,
      reference_id = record_id,
      pmid = .plant_first_pattern(flat, "\\b\\d{7,9}\\b"),
      doi = .plant_first_pattern(flat, "10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+"),
      plant_part = .plant_lotus_plant_part(flat),
      tissue = .plant_lotus_tissue(flat),
      method = .plant_lotus_method(flat),
      occurrence_type = "database_taxon_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(matched_rank == "species", "high", "medium"),
      curation_flag = "database_record_review_recommended",
      evidence_tier = ifelse(matched_rank == "species",
                             "direct_species_database",
                             paste0(matched_rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_occurrences(rows)
}

.plant_knapsack_rows = function(query_row, html, matched_taxon, matched_rank,
                                url, max_records) {
  if (is.null(html) || length(.uaf_non_empty(html)) < 1) {
    return(.plant_empty_occurrences())
  }
  table_rows = .plant_html_rows(html)
  if (length(table_rows) < 1) return(.plant_empty_occurrences())
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records)) table_rows = utils::head(table_rows, max_records)
  rows = list()
  for (cells in table_rows) {
    if (length(cells) < 2) next
    c_id = .uaf_first_non_empty_text(cells[grepl("^C\\d{6,}$", cells)])
    compound_name = .plant_knapsack_compound_name(cells, query_row)
    if (is.na(compound_name)) next
    evidence_text = .plant_truncate(paste(cells, collapse = " | "), 800)
    rows[[length(rows) + 1]] = data.frame(
      query_plant = query_row$query_plant,
      query_plant_clean = query_row$query_plant_clean,
      matched_taxon = matched_taxon,
      matched_rank = matched_rank,
      species = query_row$species,
      genus = query_row$genus,
      family = query_row$family,
      compound_name = compound_name,
      compound_name_clean = .plant_clean_compound(compound_name),
      compound_id = c_id,
      compound_id_type = ifelse(is.na(c_id), NA_character_, "KNApSAcK_C_ID"),
      source_database = "KNApSAcK",
      source_record_id = c_id,
      evidence_text = evidence_text,
      evidence_url = url,
      reference_id = c_id,
      pmid = NA_character_,
      doi = NA_character_,
      plant_part = NA_character_,
      tissue = NA_character_,
      method = NA_character_,
      occurrence_type = "organism_metabolite_record",
      retrieved_at = .plant_timestamp(),
      confidence = ifelse(matched_rank == "species", "high", "medium"),
      curation_flag = "database_record_review_recommended",
      evidence_tier = ifelse(matched_rank == "species",
                             "direct_species_database",
                             paste0(matched_rank, "_database_fallback")),
      stringsAsFactors = FALSE
    )
  }
  .plant_bind_occurrences(rows)
}

.plant_ncbi_taxonomy_search = function(species, cache, cache_dir, throttle,
                                       ncbi_email, ncbi_tool, ncbi_api_key,
                                       request_fun) {
  term = paste0('"', species, '"[Scientific Name]')
  url = .plant_ncbi_url(
    endpoint = "esearch.fcgi",
    params = c(db = "taxonomy", term = term, retmode = "json",
               retmax = "1", tool = ncbi_tool, email = ncbi_email),
    api_key = ncbi_api_key
  )
  .plant_fetch_json(url, cache, file.path(cache_dir, "ncbi_taxonomy"),
                    throttle, request_fun)
}

.plant_pubchem_taxonomy_rows = function(query_row, taxid, result, url,
                                        max_records) {
  annotations = .pubchem_parse_pugview(
    json = result,
    cid = NA_integer_,
    query = query_row$species,
    heading = paste0("Taxonomy: ", query_row$species),
    pubchem_url = url
  )
  if (!is.data.frame(annotations) || nrow(annotations) < 1) {
    return(.plant_empty_occurrences())
  }
  rows = list()
  for (i in seq_len(nrow(annotations))) {
    annotation = annotations[i, , drop = FALSE]
    ids = .plant_pubchem_compound_links(annotation)
    if (nrow(ids) < 1) next
    for (j in seq_len(nrow(ids))) {
      compound_name = .uaf_first_non_empty_text(ids$label[[j]],
                                                annotation$CleanValue,
                                                annotation$Name)
      if (is.na(compound_name) ||
          .plant_text_has_taxon(compound_name, query_row)) {
        next
      }
      rows[[length(rows) + 1]] = data.frame(
        query_plant = query_row$query_plant,
        query_plant_clean = query_row$query_plant_clean,
        matched_taxon = query_row$species,
        matched_rank = "species",
        species = query_row$species,
        genus = query_row$genus,
        family = query_row$family,
        compound_name = compound_name,
        compound_name_clean = .plant_clean_compound(compound_name),
        compound_id = ids$cid[[j]],
        compound_id_type = "PubChem_CID",
        source_database = "PubChem Taxonomy",
        source_record_id = paste("taxonomy", taxid, "cid", ids$cid[[j]],
                                 sep = ":"),
        evidence_text = .plant_truncate(annotation$CleanValue, 800),
        evidence_url = .uaf_first_non_empty_text(annotation$SourceURL,
                                                annotation$PubChemURL,
                                                url),
        reference_id = .uaf_first_non_empty_text(annotation$Source,
                                                annotation$SourceURL),
        pmid = .plant_first_pattern(annotation, "\\b\\d{7,9}\\b"),
        doi = .plant_first_pattern(annotation,
                                   "10\\.\\d{4,9}/[-._;()/:A-Za-z0-9]+"),
        plant_part = NA_character_,
        tissue = NA_character_,
        method = NA_character_,
        occurrence_type = "pubchem_taxonomy_annotation",
        retrieved_at = .plant_timestamp(),
        confidence = "medium",
        curation_flag = "source_annotation_review_recommended",
        evidence_tier = "direct_species_database",
        stringsAsFactors = FALSE
      )
    }
  }
  occurrences = .plant_bind_occurrences(rows)
  max_records = .plant_max_records(max_records)
  if (is.finite(max_records) && nrow(occurrences) > max_records) {
    occurrences = occurrences[seq_len(max_records), , drop = FALSE]
  }
  occurrences
}

.plant_article_id = function(item, id_type) {
  articleids = tryCatch(item$articleids, error = function(error) NULL)
  if (is.null(articleids)) return(NA_character_)
  for (entry in articleids) {
    if (identical(tolower(.uaf_first_non_empty_text(entry$idtype)),
                  tolower(id_type))) {
      return(.uaf_first_non_empty_text(entry$value))
    }
  }
  NA_character_
}

.plant_pubtator_rows = function(query_row, result) {
  docs = .plant_pubtator_docs(result)
  rows = list()
  for (doc in docs) {
    pmid = .uaf_first_non_empty_text(doc$pmid, doc$id, doc$sourceid)
    text = .pubchem_collapse(c(doc$title, doc$abstract,
                               doc$text_hl, doc$passages[[1]]$text))
    chemicals = unique(.uaf_non_empty(c(
      .plant_pubtator_mentions(doc, c("Chemical", "chemical")),
      .plant_pubtator_highlight_mentions(doc, "CHEMICAL")
    )))
    species_mentions = unique(.uaf_non_empty(c(
      .plant_pubtator_mentions(doc, c("Species", "species", "Organism",
                                      "Taxon")),
      .plant_pubtator_highlight_mentions(doc, "SPECIES")
    )))
    if (length(species_mentions) < 1 &&
        .plant_text_contains(text, query_row$species)) {
      species_mentions = query_row$species
    }
    if (length(species_mentions) < 1 ||
        !.plant_text_contains(.pubchem_collapse(species_mentions),
                              query_row$species)) {
      next
    }
    if (length(chemicals) < 1) chemicals = NA_character_
    for (chemical in chemicals) {
      rows[[length(rows) + 1]] = data.frame(
        query_plant = query_row$query_plant,
        query_plant_clean = query_row$query_plant_clean,
        species = query_row$species,
        genus = query_row$genus,
        family = query_row$family,
        source_database = "PubTator",
        source_record_id = pmid,
        pmid = pmid,
        doi = NA_character_,
        title = .uaf_first_non_empty_text(doc$title),
        abstract = .uaf_first_non_empty_text(doc$abstract),
        chemical_mention = chemical,
        species_mention = .pubchem_collapse(species_mentions),
        evidence_text = text,
        evidence_url = ifelse(is.na(pmid), NA_character_,
                              paste0("https://pubmed.ncbi.nlm.nih.gov/",
                                     pmid, "/")),
        retrieved_at = .plant_timestamp(),
        confidence = "low",
        curation_flag = "candidate_co_mention",
        evidence_tier = "direct_species_pubtator_candidate",
        stringsAsFactors = FALSE
      )
    }
  }
  list(LiteratureCandidates = .plant_bind_tables(rows, .plant_literature_cols()))
}

.plant_pubtator_docs = function(result) {
  if (is.null(result)) return(list())
  if (!is.null(result$PubTator3) && is.list(result$PubTator3)) {
    return(result$PubTator3)
  }
  if (!is.null(result$documents) && is.list(result$documents)) {
    return(result$documents)
  }
  if (!is.null(result$results) && is.list(result$results)) {
    return(result$results)
  }
  if (!is.null(result$passages) || !is.null(result$annotations)) {
    return(list(result))
  }
  list()
}

.plant_pubtator_mentions = function(doc, types) {
  annotations = doc$annotations
  if (is.null(annotations) && !is.null(doc$passages)) {
    annotations = unlist(lapply(doc$passages, function(passage) {
      passage$annotations
    }), recursive = FALSE)
  }
  if (is.null(annotations)) return(character())
  values = character()
  for (ann in annotations) {
    ann_type = .uaf_first_non_empty_text(ann$infons$type, ann$type,
                                         ann$infons$identifier)
    if (!ann_type %in% types) next
    values = c(values, .uaf_first_non_empty_text(ann$text))
  }
  unique(.uaf_non_empty(values))
}

.plant_pubtator_highlight_mentions = function(doc, type) {
  text = .uaf_first_non_empty_text(doc$text_hl, doc$highlight,
                                   doc$snippet)
  if (is.na(text) || text == "") return(character())
  pattern = paste0("@", type, "_[^[:space:]]+\\s+@@@(.*?)@@@")
  matches = gregexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  raw = regmatches(text, matches)[[1]]
  if (length(raw) < 1 || identical(raw, "")) return(character())
  values = vapply(raw, function(item) {
    out = sub(pattern, "\\1", item, perl = TRUE, ignore.case = TRUE)
    .plant_html_text(out)
  }, character(1))
  unique(.uaf_non_empty(values))
}

.plant_occurrences_from_pubtator = function(literature) {
  literature = .plant_normalize_literature(literature, "pubtator")
  if (nrow(literature) < 1) return(.plant_empty_occurrences())
  keep = !is.na(literature$chemical_mention) & literature$chemical_mention != ""
  literature = literature[keep, , drop = FALSE]
  if (nrow(literature) < 1) return(.plant_empty_occurrences())
  out = data.frame(
    query_plant = literature$query_plant,
    query_plant_clean = literature$query_plant_clean,
    matched_taxon = literature$species,
    matched_rank = "species",
    species = literature$species,
    genus = literature$genus,
    family = literature$family,
    compound_name = literature$chemical_mention,
    compound_name_clean = .plant_clean_compound(literature$chemical_mention),
    compound_id = NA_character_,
    compound_id_type = NA_character_,
    source_database = "PubTator",
    source_record_id = literature$source_record_id,
    evidence_text = literature$evidence_text,
    evidence_url = literature$evidence_url,
    reference_id = literature$source_record_id,
    pmid = literature$pmid,
    doi = literature$doi,
    plant_part = NA_character_,
    tissue = NA_character_,
    method = NA_character_,
    occurrence_type = "candidate_co_mention",
    retrieved_at = literature$retrieved_at,
    confidence = "low",
    curation_flag = "candidate_co_mention",
    evidence_tier = "direct_species_pubtator_candidate",
    stringsAsFactors = FALSE
  )
  .plant_normalize_occurrences(out, "pubtator")
}

.plant_json_records = function(x, depth = 0) {
  if (depth > 5 || is.null(x)) return(list())
  if (is.data.frame(x)) {
    return(lapply(seq_len(nrow(x)), function(i) {
      as.list(x[i, , drop = FALSE])
    }))
  }
  if (!is.list(x)) return(list())
  names_x = .plant_normalize_column_names(names(x))
  record_keys = c("name", "title", "compound", "molecule", "smiles", "inchi",
                  "inchikey", "lotus", "wikidata", "organism", "taxon",
                  "species", "genus", "family")
  if (any(grepl(paste(record_keys, collapse = "|"), names_x,
                ignore.case = TRUE))) {
    return(list(x))
  }
  out = list()
  for (item in x) {
    out = c(out, .plant_json_records(item, depth + 1))
  }
  out
}

.plant_flatten_record = function(x, prefix = character(), depth = 0) {
  if (depth > 5 || is.null(x)) return(character())
  if (is.data.frame(x)) {
    if (nrow(x) < 1) return(character())
    x = as.list(x[1, , drop = FALSE])
  }
  if (is.list(x)) {
    out = character()
    names_x = names(x)
    if (is.null(names_x)) names_x = paste0("item", seq_along(x))
    for (i in seq_along(x)) {
      name = .plant_normalize_column_names(names_x[[i]])
      if (is.na(name) || name == "") name = paste0("item", i)
      child = .plant_flatten_record(x[[i]], c(prefix, name), depth + 1)
      out = c(out, child)
    }
    return(out)
  }
  value = .uaf_squish_text(as.character(unlist(x, use.names = FALSE)))
  value = .uaf_non_empty(value)
  if (length(value) < 1) return(character())
  key = paste(prefix, collapse = "_")
  if (key == "") key = "value"
  stats::setNames(value, rep(key, length(value)))
}

.plant_lotus_compound_name = function(flat, query_row) {
  priority = c("compound_name", "molecule_name", "natural_product_name",
               "traditional_name", "iupac_name", "name", "title")
  keys = names(flat)
  for (pattern in priority) {
    values = flat[grepl(pattern, keys, ignore.case = TRUE)]
    values = .uaf_non_empty(values)
    values = values[!vapply(values, .plant_text_has_taxon, logical(1),
                            query_row = query_row)]
    values = values[!grepl("^LTS\\d+|^Q\\d+$", values, ignore.case = TRUE)]
    if (length(values) > 0) return(values[[1]])
  }
  NA_character_
}

.plant_lotus_record_id = function(flat) {
  values = flat[grepl("lotus|wikidata|record|identifier|^id$",
                      names(flat), ignore.case = TRUE)]
  values = .uaf_non_empty(values)
  hit = values[grepl("^LTS\\d+|^Q\\d+$|^https?://", values,
                     ignore.case = TRUE)]
  .uaf_first_non_empty_text(hit, values)
}

.plant_knapsack_compound_name = function(cells, query_row) {
  cells = .uaf_non_empty(.uaf_squish_text(cells))
  if (length(cells) < 1) return(NA_character_)
  reject = grepl(
    paste(c("^C\\d{6,}$", "^CAS$", "^C_ID$", "metabolite", "organism",
            "formula", "molecular", "^mw$", "exact", "inchi", "smiles",
            "search", "image", "link", "name", "taxonomy", "all rights"),
          collapse = "|"),
    cells,
    ignore.case = TRUE
  )
  taxon_cell = vapply(cells, .plant_text_has_taxon, logical(1),
                      query_row = query_row)
  reject = reject | .plant_formula_like(cells) | taxon_cell
  candidates = cells[!reject]
  candidates = candidates[grepl("[A-Za-z]", candidates)]
  candidates = candidates[nchar(candidates) >= 3]
  .uaf_first_non_empty_text(candidates)
}

.plant_html_rows = function(html) {
  html = paste(as.character(html), collapse = "\n")
  rows = regmatches(html, gregexpr("<tr[^>]*>.*?</tr>", html,
                                   perl = TRUE, ignore.case = TRUE))[[1]]
  out = list()
  if (length(rows) > 0 && !identical(rows, "")) {
    for (row in rows) {
      cells = regmatches(row, gregexpr("<t[dh][^>]*>.*?</t[dh]>", row,
                                       perl = TRUE, ignore.case = TRUE))[[1]]
      cells = .uaf_non_empty(.plant_html_text(cells))
      if (length(cells) > 0) out[[length(out) + 1]] = cells
    }
  }
  if (length(out) > 0) return(out)
  text_rows = strsplit(.plant_html_text(html), "\n", fixed = TRUE)[[1]]
  text_rows = .uaf_non_empty(text_rows)
  for (row in text_rows) {
    cells = .uaf_non_empty(strsplit(row, "\\s{2,}|\\t", perl = TRUE)[[1]])
    if (length(cells) > 1) out[[length(out) + 1]] = cells
  }
  out
}

.plant_html_text = function(x) {
  x = gsub("<script[^>]*>.*?</script>", " ", x, ignore.case = TRUE,
           perl = TRUE)
  x = gsub("<style[^>]*>.*?</style>", " ", x, ignore.case = TRUE,
           perl = TRUE)
  x = gsub("<[^>]+>", " ", x, perl = TRUE)
  x = gsub("&nbsp;", " ", x, ignore.case = TRUE)
  x = gsub("&amp;", "&", x, ignore.case = TRUE)
  x = gsub("&quot;", "\"", x, ignore.case = TRUE)
  x = gsub("&#39;|&apos;", "'", x, ignore.case = TRUE)
  .uaf_squish_text(x)
}

.plant_pubchem_compound_links = function(annotation) {
  urls = .pubchem_split_collapsed(annotation$MarkupURL)
  labels = .pubchem_split_collapsed(annotation$MarkupText)
  rows = list()
  for (i in seq_along(urls)) {
    cid = .plant_pubchem_cid_from_url(urls[[i]])
    if (is.na(cid)) next
    label = if (length(labels) >= i) labels[[i]] else NA_character_
    rows[[length(rows) + 1]] = data.frame(
      cid = cid,
      label = label,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) < 1) return(.uaf_empty_table(c("cid", "label")))
  unique(do.call(rbind, rows))
}

.plant_pubchem_cid_from_url = function(url) {
  url = .uaf_first_non_empty_text(url)
  if (is.na(url)) return(NA_character_)
  patterns = c("/compound/(\\d+)", "[?&]cid=(\\d+)", "/cid/(\\d+)",
               "CID(\\d+)")
  for (pattern in patterns) {
    hit = regexec(pattern, url, ignore.case = TRUE)
    match = regmatches(url, hit)[[1]]
    if (length(match) > 1) return(match[[2]])
  }
  NA_character_
}

.plant_matched_rank_from_text = function(query_row, text) {
  if (.plant_text_contains(text, query_row$species)) return("species")
  if (.plant_text_contains(text, query_row$genus)) return("genus")
  if (.plant_text_contains(text, query_row$family)) return("family")
  NA_character_
}

.plant_text_has_taxon = function(text, query_row) {
  .plant_text_contains(text, query_row$species) ||
    .plant_text_contains(text, query_row$genus) ||
    .plant_text_contains(text, query_row$family)
}

.plant_text_contains = function(text, value) {
  text = .plant_clean_name(.pubchem_collapse(text))
  value = .plant_clean_name(value)
  value = .uaf_non_empty(value)
  if (length(value) < 1 || is.na(text) || text == "") return(FALSE)
  any(vapply(value, function(v) {
    !is.na(v) && v != "" && grepl(v, text, fixed = TRUE)
  }, logical(1)))
}

.plant_formula_like = function(x) {
  x = .uaf_squish_text(x)
  grepl("^[A-Z][a-z]?\\d*([A-Z][a-z]?\\d*)+$", x) |
    grepl("^[0-9.]+$", x)
}

.plant_first_named = function(flat, pattern) {
  values = flat[grepl(pattern, names(flat), ignore.case = TRUE)]
  .uaf_first_non_empty_text(values)
}

.plant_first_named_filtered = function(flat, include, exclude = NA_character_) {
  nm = names(flat)
  keep = grepl(include, nm, ignore.case = TRUE, perl = TRUE)
  if (!is.na(exclude)) {
    keep = keep & !grepl(exclude, nm, ignore.case = TRUE, perl = TRUE)
  }
  values = flat[keep]
  if (length(values) < 1) return(NA_character_)
  for (value in unname(values)) {
    value = .uaf_first_non_empty_text(value)
    if (is.na(value) || value == "") next
    if (grepl("^[0-9]+$", value)) next
    if (grepl("^https?://", value, ignore.case = TRUE)) next
    return(value)
  }
  NA_character_
}

.plant_lotus_context_exclude = function() {
  paste(c("organism", "taxonomy", "taxon", "kingdom", "phylum", "classx",
          "family", "genus", "species", "cleaned_organism_id", "wikidata",
          "reference", "pubchemfingerprint", "fragment",
          "chemicaltaxonomy"), collapse = "|")
}

.plant_lotus_plant_part = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])plant[_ ]?part([._]|$)",
            "(^|[._])plantpart([._]|$)",
            "(^|[._])part[_ ]?used([._]|$)",
            "(^|[._])partused([._]|$)",
            "(^|[._])sample[_ ]?organ([._]|$)",
            "(^|[._])sample[_ ]?part([._]|$)",
            "(^|[._])material([._]|$)"), collapse = "|"),
    .plant_lotus_context_exclude()
  )
}

.plant_lotus_tissue = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])tissue([._]|$)",
            "(^|[._])sample[_ ]?tissue([._]|$)"), collapse = "|"),
    .plant_lotus_context_exclude()
  )
}

.plant_lotus_method = function(flat) {
  .plant_first_named_filtered(
    flat,
    paste(c("(^|[._])method([._]|$)",
            "(^|[._])analysis([._]|$)",
            "(^|[._])technique([._]|$)",
            "(^|[._])instrument([._]|$)"), collapse = "|"),
    "chemicaltaxonomy|pubchemfingerprint|fragment"
  )
}

.plant_first_pattern = function(values, pattern) {
  text = .pubchem_collapse(unlist(values, use.names = FALSE))
  hit = regexpr(pattern, text, perl = TRUE, ignore.case = TRUE)
  if (hit[[1]] < 0) return(NA_character_)
  regmatches(text, hit)
}

.plant_max_records = function(max_records) {
  max_records = suppressWarnings(as.numeric(max_records[[1]]))
  if (!is.finite(max_records) || max_records < 1) return(Inf)
  max_records
}

.plant_truncate = function(x, n) {
  x = .uaf_first_non_empty_text(x)
  if (is.na(x) || nchar(x) <= n) return(x)
  paste0(substr(x, 1, n - 3), "...")
}

.plant_normalize_column_names = function(x) {
  x = gsub("[^A-Za-z0-9]+", "_", x)
  x = gsub("_+", "_", x)
  x = gsub("^_|_$", "", x)
  tolower(x)
}

.plant_col_or_default = function(table, column, default) {
  n = if (is.data.frame(table)) nrow(table) else 0
  if (!is.data.frame(table) || !column %in% names(table)) {
    return(rep(default, n))
  }
  value = table[[column]]
  value = .uaf_squish_text(value)
  value[is.na(value) | value == ""] = default
  value
}

.plant_clean_name = function(x) {
  x = tolower(.uaf_squish_text(x))
  x = gsub("\\s+", " ", x)
  x
}

.plant_clean_compound = function(x) {
  x = tolower(.uaf_squish_text(x))
  x = gsub("[^a-z0-9]+", "_", x)
  x = gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

.plant_slug = function(x) {
  x = .plant_clean_compound(x)
  x[is.na(x) | x == ""] = NA_character_
  x
}

.plant_genus = function(species) {
  species = .uaf_squish_text(species)
  vapply(strsplit(species, "\\s+"), function(parts) {
    if (length(parts) < 1) return(NA_character_)
    parts[[1]]
  }, character(1))
}

.plant_is_binomial = function(x) {
  grepl("^[A-Z][A-Za-z-]+\\s+[a-z][A-Za-z-]+", .uaf_squish_text(x))
}

.plant_cache_dir = function(cache_dir) {
  if (!is.null(cache_dir) && length(.uaf_non_empty(cache_dir)) > 0) {
    return(.uaf_non_empty(cache_dir)[[1]])
  }
  out = tryCatch(tools::R_user_dir("uafR", "cache"),
                 error = function(error) NULL)
  if (is.null(out)) out = file.path(tempdir(), "uafR-plant-cache")
  file.path(out, "plant_phytochemistry")
}

.plant_timestamp = function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

.plant_normalize_evidence_tier = function(evidence_tier, source_database,
                                          matched_rank) {
  evidence_tier = .uaf_squish_text(evidence_tier)
  source_database = tolower(.uaf_squish_text(source_database))
  matched_rank = tolower(.uaf_squish_text(matched_rank))
  out = evidence_tier
  missing = is.na(out) | out == ""
  out[missing & source_database %in% c("manual", "curated")] = "manual_curated"
  out[missing & source_database %in% c("pubtator")] =
    "direct_species_pubtator_candidate"
  out[missing & source_database %in% c("pubmed")] = "direct_species_literature"
  out[missing & matched_rank == "genus"] = "genus_database_fallback"
  out[missing & matched_rank == "family"] = "family_database_fallback"
  out[missing] = "direct_species_database"
  out[!out %in% .plant_evidence_tiers()] = "unresolved"
  out
}

.plant_evidence_tier_confidence = function(evidence_tier) {
  tier = .plant_normalize_evidence_tier(evidence_tier, "", "species")
  ifelse(tier %in% c("direct_species_database", "manual_curated"), "high",
         ifelse(tier %in% c("direct_species_literature",
                            "genus_database_fallback",
                            "genus_literature_fallback"), "medium", "low"))
}

.plant_normalize_confidence = function(confidence, evidence_tier) {
  confidence = tolower(.uaf_squish_text(confidence))
  fallback = .plant_evidence_tier_confidence(evidence_tier)
  confidence[is.na(confidence) | confidence == ""] =
    fallback[is.na(confidence) | confidence == ""]
  confidence[!confidence %in% .plant_confidence_values()] = "unknown"
  confidence
}

.plant_confidence_score = function(confidence) {
  confidence = tolower(.uaf_squish_text(confidence))
  if (length(confidence) == 1 && confidence %in% c("low", "medium", "high",
                                                   "unknown")) {
    return(switch(confidence, low = 0.35, medium = 0.65, high = 1,
                  unknown = 0))
  }
  out = suppressWarnings(as.numeric(confidence))
  out[confidence == "low"] = 0.35
  out[confidence == "medium"] = 0.65
  out[confidence == "high"] = 1
  out[confidence == "unknown" | is.na(confidence)] = 0
  out[is.na(out)] = 0
  out
}

.plant_clean_context_value = function(x) {
  out = .uaf_squish_text(x)
  lower = tolower(out)
  out[is.na(out) | out == "" |
        lower %in% c("na", "n/a", "none", "unknown", "unspecified",
                     "not available", "not reported")] = NA_character_
  out[grepl("^[0-9]+$", out)] = NA_character_
  out[grepl("^https?://", out, ignore.case = TRUE)] = NA_character_
  out
}

.plant_context_group = function(value, supplied = NA_character_,
                                kind = c("plant_part", "tissue")) {
  kind = match.arg(kind)
  supplied = .plant_matrix_key(supplied)
  valid = .plant_context_group_values()
  supplied[supplied == "unknown"] = NA_character_
  supplied[!supplied %in% valid] = NA_character_
  computed = vapply(.plant_clean_context_value(value), function(item) {
    text = tolower(.uaf_first_non_empty_text(item))
    if (is.na(text) || text == "") return("unknown")
    if (grepl("root|rhizome|tuber|bulb|belowground", text)) {
      return("root_belowground")
    }
    if (grepl("rhizosphere|exudate|exudation", text)) {
      return("exudate_rhizosphere")
    }
    if (grepl("leaf|leaves|foliar|needle", text)) return("leaf")
    if (grepl("stem|shoot|twig|branch|culm", text)) return("stem_shoot")
    if (grepl("bark|wood|xylem|phloem", text)) return("bark_wood")
    if (grepl("flower|floral|inflorescence|petal|pollen|stamen", text)) {
      return("flower")
    }
    if (grepl("seedling|seedlings", text)) return("whole_plant")
    if (grepl("fruit|\\bseed\\b|\\bseeds\\b|grain|kernel|berry|pod|achene|nut", text)) {
      return("fruit_seed")
    }
    if (grepl("aerial|aboveground|above-ground", text)) return("aerial")
    if (grepl("whole plant|whole-plant|whole organism|entire plant", text)) {
      return("whole_plant")
    }
    if (grepl("vascular|vein|vasculature", text)) return("vascular")
    if (grepl("gland|secretory|trichome|resin duct", text)) return("secretory")
    if (grepl("epiderm|cuticle|surface", text)) return("epidermal")
    if (grepl("microb|endophy|rhizob", text)) return("microbial")
    if (grepl("extract|essential oil|oil", text)) return("extract_unspecified")
    "other"
  }, character(1))
  out = ifelse(!is.na(supplied) & supplied != "", supplied, computed)
  out[!out %in% valid] = "unknown"
  out
}

.plant_method_group = function(value, supplied = NA_character_,
                               source_database = NA_character_,
                               evidence_tier = NA_character_) {
  supplied = .plant_matrix_key(supplied)
  valid = .plant_method_group_values()
  supplied[supplied == "unknown"] = NA_character_
  supplied[!supplied %in% valid] = NA_character_
  computed = vapply(.plant_clean_context_value(value), function(item) {
    text = tolower(.uaf_first_non_empty_text(item))
    if (is.na(text) || text == "") return("unknown")
    if (grepl("gc[- ]?ms|gas chromat", text)) return("gc_ms")
    if (grepl("lc[- ]?ms|liquid chromat.*mass", text)) return("lc_ms")
    if (grepl("hplc|uplc", text)) return("hplc")
    if (grepl("\\bnmr\\b|nuclear magnetic resonance", text)) return("nmr")
    if (grepl("mass spectrom|\\bms\\b", text)) return("mass_spectrometry")
    if (grepl("chromat", text)) return("chromatography")
    if (grepl("spectroscop|uv-vis|infrared|ftir", text)) {
      return("spectroscopy")
    }
    if (grepl("literature|paper|publication|curat", text)) {
      return("literature_curation")
    }
    "other"
  }, character(1))
  missing_method = computed == "unknown"
  source = tolower(.uaf_squish_text(source_database))
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        "species")
  computed[missing_method & source %in% c("lotus", "knapsack", "npass",
                                          "pubchem taxonomy")] =
    "database_record"
  computed[missing_method & grepl("literature|pubtator", tier)] =
    "literature_curation"
  out = ifelse(!is.na(supplied) & supplied != "", supplied, computed)
  out[!out %in% valid] = "unknown"
  out
}

.plant_context_status = function(plant_part_group, tissue_group, method_group) {
  part_known = !.plant_matrix_key(plant_part_group) %in%
    c("unknown", "extract_unspecified")
  tissue_known = !.plant_matrix_key(tissue_group) %in%
    c("unknown", "extract_unspecified")
  method_known = !.plant_matrix_key(method_group) %in%
    c("unknown", "database_record", "literature_curation")
  ifelse((part_known | tissue_known) & method_known,
         "plant_part_and_method_known",
         ifelse(part_known | tissue_known, "plant_part_known",
                ifelse(method_known, "method_known", "context_missing")))
}

.plant_occurrence_status = function(evidence_tier, matched_rank,
                                    source_database, curation_flag) {
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        matched_rank)
  rank = tolower(.uaf_squish_text(matched_rank))
  status = rep("unknown", length(tier))
  status[tier == "manual_curated"] = "curated_reported"
  status[tier == "direct_species_database" & rank == "species"] =
    "direct_reported"
  status[tier == "direct_species_literature"] = "literature_reported"
  status[tier == "direct_species_pubtator_candidate"] = "candidate"
  status[grepl("fallback", tier)] = "taxon_fallback"
  status[tier == "unresolved"] = "unresolved"
  status
}

.plant_occurrence_basis = function(evidence_tier, matched_rank,
                                   source_database) {
  tier = .plant_normalize_evidence_tier(evidence_tier, source_database,
                                        matched_rank)
  rank = tolower(.uaf_squish_text(matched_rank))
  basis = rep("unknown", length(tier))
  basis[tier == "manual_curated"] = "species_curated_record"
  basis[tier == "direct_species_database" & rank == "species"] =
    "species_database_record"
  basis[tier == "direct_species_literature"] = "species_literature_record"
  basis[tier == "direct_species_pubtator_candidate"] = "candidate_co_mention"
  basis[tier == "genus_database_fallback" | rank == "genus"] =
    "genus_fallback_record"
  basis[tier == "family_database_fallback" | rank == "family"] =
    "family_fallback_record"
  basis[tier == "unresolved"] = "unresolved"
  basis
}

.plant_evidence_quality_score = function(occurrence_status, occurrence_basis,
                                         confidence,
                                         biological_context_status) {
  status = .plant_matrix_key(occurrence_status)
  base = rep(0, length(status))
  base[status == "curated_reported"] = 0.95
  base[status == "direct_reported"] = 0.85
  base[status == "literature_reported"] = 0.55
  base[status == "taxon_fallback"] = 0.45
  base[status == "candidate"] = 0.20
  base[status == "unresolved"] = 0
  conf = .plant_confidence_score(confidence)
  context = .plant_matrix_key(biological_context_status)
  bonus = rep(0, length(status))
  bonus[context == "plant_part_and_method_known"] = 0.05
  bonus[context == "plant_part_known"] = 0.03
  bonus[context == "method_known"] = 0.02
  round(pmin(1, base * conf + bonus), 3)
}

.plant_analysis_ready = function(occurrence_status, confidence,
                                 evidence_quality_score) {
  status = .plant_matrix_key(occurrence_status)
  conf_ok = .plant_confidence_score(confidence) >= .plant_confidence_score("medium")
  score = suppressWarnings(as.numeric(evidence_quality_score))
  ifelse(status %in% c("direct_reported", "curated_reported") &
           conf_ok & !is.na(score) & score >= 0.50, "Yes", "No")
}

.plant_matrix_key = function(x) {
  key = .plant_clean_compound(x)
  key[is.na(key) | key == ""] = "unknown"
  key
}

.plant_count_summary = function(x) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(NA_character_)
  tab = sort(table(x), decreasing = TRUE)
  paste(paste(names(tab), as.integer(tab), sep = "="), collapse = "; ")
}

.plant_top_terms = function(x, n = 5) {
  x = .uaf_non_empty(x)
  if (length(x) < 1) return(NA_character_)
  tab = sort(table(x), decreasing = TRUE)
  .pubchem_collapse(utils::head(names(tab), n))
}

.plant_source_coverage_score = function(source_count, provider_diagnostics) {
  denom = if (is.data.frame(provider_diagnostics) &&
              nrow(provider_diagnostics) > 0) {
    max(1, sum(provider_diagnostics$enabled == "Yes"))
  } else {
    max(1, source_count)
  }
  round(min(1, source_count / denom), 3)
}

.plant_categorate_validation_status = function(categorate_result) {
  if (is.list(categorate_result) &&
      is.data.frame(categorate_result$ValidationSummary) &&
      "Status" %in% names(categorate_result$ValidationSummary)) {
    return(.uaf_first_non_empty_text(categorate_result$ValidationSummary$Status))
  }
  if (is.null(categorate_result)) return("not_run")
  "not_available"
}

"%||%" = function(x, y) if (is.null(x)) y else x
