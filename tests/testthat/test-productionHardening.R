test_that("workflow guide and API stability expose production classifications", {
  guide = uafRWorkflowGuide()
  stability = uafRApiStability()

  expect_true(all(c("user_has", "recommended_workflow", "primary_functions",
                    "stability") %in% names(guide)))
  expect_true("runPlantChemistryProject" %in% stability$function_name)
  expect_true(all(stability$stability %in%
                    c("stable", "experimental", "internal",
                      "project_specific")))
})

test_that("schema metadata and provider diagnostics are standardized", {
  meta = uafRSchemaMetadata("unit_test", list(alpha = 1))
  expect_true(all(c("uafR_schema_version", "uafR_package_version",
                    "created_at", "workflow_name",
                    "workflow_parameters") %in% names(meta)))
  expect_equal(meta$workflow_name, "unit_test")
  expect_equal(meta$uafR_schema_version, "1.0.0")

  diag = standardizeProviderDiagnostics(data.frame(
    provider = "PubChem",
    enabled = "Yes",
    queried = "Yes",
    request_count = 3,
    cache_hit_count = 1,
    record_count = 0,
    error_count = 0,
    message = "no records",
    stringsAsFactors = FALSE
  ))

  expect_equal(diag$query_count, 3)
  expect_equal(diag$cache_hit_count, 1)
  expect_equal(diag$no_hit_reason,
               "no_records_found_or_source_not_applicable")
})

test_that("provider contracts and claim guidance document conservative use", {
  contracts = uafRProviderContracts(c("PubChem", "PubTator"))
  claims = uafRClaimGuidance(c("tanimoto_similarity",
                               "pubtator_pubmed_candidate"))

  expect_equal(sort(contracts$provider), c("PubChem", "PubTator"))
  expect_match(contracts$no_hit_interpretation[[1]], "not evidence")
  expect_true(all(c("source_url", "fields_extracted", "terms_review_status",
                    "redistribution_default", "contract_reviewed_on") %in%
                    names(contracts)))
  expect_true(all(grepl("^https://", contracts$source_url)))
  expect_true(all(contracts$terms_review_status ==
                    "not_asserted_check_current_provider_terms"))
  expect_true(all(c("can_support", "cannot_support",
                    "recommended_safe_wording") %in% names(claims)))
  expect_true("pubtator_pubmed_candidate_only" %in%
                plantOccurrenceEvidenceDictionary()$evidence_grade)
})

test_that("cache inspection and run planning work without network", {
  cache_dir = tempfile("uafr_cache_")
  dir.create(file.path(cache_dir, "pubchem"), recursive = TRUE)
  writeLines("{}", file.path(cache_dir, "pubchem", "one.json"))
  saveRDS(list(ok = TRUE), file.path(cache_dir, "categorate_batch_0001.rds"))

  inventory = inspectUafRCache(cache_dir)
  summary = summarizeUafRCache(cache_dir)
  lotus_index = tempfile("lotus_plan_index_")
  dir.create(lotus_index)
  plan = planPlantChemistryRun(
    plants = paste("Plant", seq_len(100)),
    compounds = data.frame(compound_id = paste0("c", seq_len(25))),
    sources = c("lotus", "pubmed"),
    cache_dir = cache_dir,
    lotus_index = lotus_index,
    write_full_pairwise = TRUE
  )

  expect_equal(nrow(inventory), 2)
  expect_true("PubChem" %in% summary$provider)
  expect_equal(plan$Summary$species_count, 100)
  expect_equal(plan$OutputEstimates$estimated_unique_compounds, 25)
  expect_true(any(plan$ProviderPlan$recommended_mode == "local_index"))
  expect_true(all(c("ReadinessChecks", "RunConfiguration") %in% names(plan)))
})

test_that("large-run planning exposes a staged 705-species operating plan", {
  plants = paste("Simulata species", sprintf("%04d", seq_len(705)))
  cache_dir = tempfile("uafr_705_cache_")
  lotus_index = tempfile("uafr_705_lotus_")
  dir.create(cache_dir)
  dir.create(lotus_index)

  plan = planPlantChemistryRun(
    plants = plants,
    sources = c("lotus", "pubmed"),
    cache_dir = cache_dir,
    lotus_index = lotus_index,
    species_chunk_size = 25,
    compound_batch_size = 25,
    max_pubmed_records = 10
  )

  expect_equal(plan$Summary$species_count, 705)
  expect_equal(plan$Summary$planned_discovery_chunk_count, 29)
  expect_equal(plan$OutputEstimates$estimated_species_pair_rows,
               choose(705, 2))
  expect_equal(plan$OutputEstimates$estimated_network_request_lower_bound,
               705 * 2)
  expect_equal(plan$Summary$readiness_status, "staged_run_required")
  expect_false(any(plan$ReadinessChecks$status == "fail"))
  expect_true(any(plan$RunConfiguration$stage == "3_identity_resolution"))

  unsafe = planPlantChemistryRun(
    plants = plants,
    sources = "lotus",
    species_chunk_size = 100,
    compound_batch_size = 100,
    write_full_pairwise = TRUE
  )
  expect_equal(unsafe$Summary$readiness_status, "not_ready")
  expect_true(any(unsafe$ReadinessChecks$status == "fail"))
})

test_that("large-run planning audits species-level names before launch", {
  cache_dir = tempfile("uafr_name_audit_cache_")
  lotus_index = tempfile("uafr_name_audit_lotus_")
  dir.create(cache_dir)
  dir.create(lotus_index)
  plants = c("Acer rubrum", "Acer", "Artocarpus sp.",
             "Acer rubrum", "")

  plan = planPlantChemistryRun(
    plants = plants,
    sources = "lotus",
    cache_dir = cache_dir,
    lotus_index = lotus_index
  )

  expect_equal(plan$Summary$input_name_count, 5)
  expect_equal(plan$Summary$species_count, 3)
  expect_equal(plan$Summary$parsed_species_count, 1)
  expect_equal(plan$Summary$review_required_name_count, 2)
  expect_equal(plan$Summary$blank_input_count, 1)
  expect_equal(plan$Summary$duplicate_input_count, 1)
  expect_equal(nrow(plan$InputNameAudit), 5)
  expect_equal(plan$Summary$readiness_status, "not_ready")
  expect_equal(
    plan$ReadinessChecks$status[
      plan$ReadinessChecks$check == "species_level_name_resolution"
    ],
    "fail"
  )
  expect_setequal(
    plan$PlantQueries$species[plan$PlantQueries$review_required],
    c("Acer", "Artocarpus sp.")
  )
})

test_that("chemistry classification overrides are reproducible", {
  dictionary = chemistryComparisonDictionary()
  expect_true(all(c("volatile_specialized_metabolites",
                    "plant_hormone_signaling",
                    "xenobiotic_or_contaminant") %in%
                  dictionary$comparison_scope))
  expect_false(any(dictionary$comparison_scope %in%
                     c("volatile_specialized", "hormones_signaling")))
  expect_true(all(dictionary$metabolism_domain %in%
                    .plant_metabolism_domain_values()))

  comparability = data.frame(
    compound_id = c("c1", "c2"),
    compound_name = c("alpha", "beta"),
    compound_name_clean = c("alpha", "beta"),
    comparison_scope = c("unknown", "unknown"),
    comparison_group = c("unknown", "unknown"),
    comparable_for_matrix = c("No", "No"),
    stringsAsFactors = FALSE
  )
  overrides = data.frame(
    compound_id = "c1",
    comparison_scope = "specialized_metabolites",
    comparison_group = "terpenoids",
    comparison_subgroup = "monoterpenoids",
    metabolism_domain = "specialized",
    biosynthetic_family = "terpenoid",
    chemical_behavior = "volatile_or_semivolatile",
    stringsAsFactors = FALSE
  )

  out = applyChemistryClassificationOverrides(comparability, overrides)

  expect_equal(out$comparison_scope[[1]], "specialized_metabolites")
  expect_equal(out$comparison_group[[1]], "terpenoid")
  expect_equal(out$metabolism_domain[[1]], "specialized_metabolism")
  expect_equal(out$chemical_behavior[[1]], "volatile_semivolatile")
  expect_equal(out$classification_source[[1]], "user_override")
  expect_equal(out$comparable_for_matrix[[1]], "Yes")
  expect_equal(out$comparison_scope[[2]], "unknown")

  expect_error(
    standardizeChemistryClassificationOverrides(data.frame(
      compound_id = "c2",
      comparison_scope = "primary_metabolites",
      comparison_group = "volatile_terpenoid",
      stringsAsFactors = FALSE
    )),
    "Inconsistent comparison scope/group"
  )
  expect_error(
    standardizeChemistryClassificationOverrides(data.frame(
      compound_id = "c2",
      comparison_scope = "invented_scope",
      stringsAsFactors = FALSE
    )),
    "Invalid `comparison_scope`"
  )
})

test_that("evidence filters and expanded feature matrices are analysis-ready", {
  membership = data.frame(
    species = c("Plant A", "Plant A", "Plant B"),
    compound_id = c("c1", "c2", "c3"),
    compound_name = c("alpha", "beta", "gamma"),
    source_database = c("LOTUS", "PubTator", "PubMed"),
    source_record_id = c("L1", "", "P1"),
    evidence_tier = c("direct_species_database",
                      "direct_species_pubtator_candidate",
                      "direct_species_literature"),
    matched_rank = "species",
    SMILES = c("CCO", NA, "CCC"),
    comparison_scope = c("specialized_metabolites", "unknown",
                         "primary_metabolites"),
    comparison_group = c("terpenoids", "unknown", "amino_acids"),
    comparable_for_matrix = c("Yes", "No", "Yes"),
    plant_part_group = c("leaf", "unknown", "root"),
    tissue_group = c("leaf", "unknown", "root"),
    method_group = c("LC-MS", "unknown", "GC-MS"),
    confidence = c("high", "low", "medium"),
    stringsAsFactors = FALSE
  )

  direct = filterPlantEvidenceDirect(membership)
  comparable = filterPlantEvidenceComparable(membership)
  review = filterPlantEvidenceReviewRequired(membership)
  features = exportPlantChemistryFeatureSet(
    membership,
    modes = c("count", "binary", "fraction", "confidence"),
    species_universe = c("Plant A", "Plant B", "Plant C")
  )

  expect_equal(nrow(direct), 2)
  expect_equal(nrow(comparable), 2)
  expect_equal(nrow(review), 1)
  expect_true(all(c("EvidenceGradeCountMatrix",
                    "PlantPartCountMatrix",
                    "ComparisonGroupConfidenceMatrix") %in% names(features)))
  expect_true("confidence" %in% features$Manifest$Mode)
  feature_tables = features$Manifest$Table
  expect_true(all(vapply(feature_tables, function(table) {
    identical(features[[table]]$species_id, features$SpeciesMetadata$species_id)
  }, logical(1))))
  expect_equal(nrow(features$SpeciesMetadata), 3)
  expect_equal(features$SpeciesMetadata$chemistry_record_status,
               c("records_present", "records_present",
                 "no_records_in_membership"))
  plant_c = features$ComparisonGroupCountMatrix$species == "Plant C"
  expect_equal(sum(features$ComparisonGroupCountMatrix[plant_c, -c(1, 2),
                                                        drop = FALSE]), 0)
})

test_that("finalized bundle manifests include schema metadata", {
  out_dir = tempfile("uafr_schema_bundle_")
  plant_compounds = data.frame(
    species = c("Plant A", "Plant B"),
    compound_name = c("alpha", "beta"),
    compound_id = c("c1", "c2"),
    source_database = "fixture",
    source_record_id = c("F1", "F2"),
    evidence_tier = "direct_species_database",
    matched_rank = "species",
    SMILES = c("CCO", "CCC"),
    comparison_scope = "specialized_metabolites",
    comparison_group = "terpenoids",
    comparable_for_matrix = "Yes",
    stringsAsFactors = FALSE
  )

  manifest = runPlantChemistryProject(
    plant_list = c("Plant A", "Plant B"),
    output_dir = out_dir,
    plant_compounds = plant_compounds,
    project_id = "schema_test",
    overwrite = TRUE
  )
  bundle_manifest = utils::read.csv(
    file.path(manifest$Project$bundle_dir[[1]], "01_ExportManifest.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  expect_true(all(c("uafR_schema_version", "uafR_package_version",
                    "workflow_name", "workflow_parameters", "OutputPath",
                    "PathType", "ArtifactPath", "FileSizeBytes",
                    "ChecksumAlgorithm", "ArtifactChecksum",
                    "ChecksumStatus") %in%
                    names(bundle_manifest)))
  expect_true(all(bundle_manifest$uafR_schema_version == "1.0.0"))
  expect_true(all(bundle_manifest$OutputPath == "."))
  expect_true(all(bundle_manifest$PathType == "bundle_relative"))
  expect_true(all(bundle_manifest$ArtifactPath == bundle_manifest$FileName))
  expect_true(all(bundle_manifest$ChecksumStatus[
    bundle_manifest$Table != "ExportManifest"
  ] == "computed"))
  expect_true("FeatureMatrixManifest" %in% bundle_manifest$Table)
  validation = validatePlantChemistryAnalysisBundle(
    manifest$Project$bundle_dir[[1]]
  )
  expect_equal(validation$Summary$ExportReadyStatus, "pass")
  expect_equal(validation$Summary$ArtifactWarnCount, 0)
  expect_equal(validation$Summary$ArtifactFailCount, 0)
  expect_true(all(validation$ArtifactValidation$Status == "pass"))
  expect_true(all(validation$ArtifactValidation$PortablePath == "Yes"))
  expect_equal(validation$Summary$ManifestReferenceFailCount, 0)
  expect_true(all(validation$ManifestReferences$Status == "pass"))
  expect_true(all(validation$ManifestReferences$FileExists == "Yes"))

  artifact_row = which(bundle_manifest$Table ==
                         "PlantCompoundMembershipEnriched")[[1]]
  bundle_manifest$ChecksumAlgorithm[[artifact_row]] = "SHA256"
  .categorate_write_csv_file(
    bundle_manifest,
    file.path(manifest$Project$bundle_dir[[1]], "01_ExportManifest.csv")
  )
  invalid_algorithm = validatePlantChemistryAnalysisBundle(
    manifest$Project$bundle_dir[[1]]
  )
  expect_equal(invalid_algorithm$Summary$ExportReadyStatus, "fail")
  expect_match(invalid_algorithm$ArtifactValidation$Message[
    invalid_algorithm$ArtifactValidation$Table ==
      "PlantCompoundMembershipEnriched"
  ], "algorithm")
})

test_that("bundle validation fails broken nested feature references", {
  out_dir = tempfile("uafr_broken_feature_manifest_")
  manifest = runPlantChemistryProject(
    plant_list = c("Plant A", "Plant B"),
    output_dir = out_dir,
    plant_compounds = data.frame(
      species = c("Plant A", "Plant B"),
      compound_name = c("alpha", "beta"),
      compound_id = c("c1", "c2"),
      source_database = "fixture",
      source_record_id = c("F1", "F2"),
      evidence_tier = "direct_species_database",
      matched_rank = "species",
      SMILES = c("CCO", "CCC"),
      comparison_scope = "specialized_metabolites",
      comparison_group = "terpenoid",
      comparable_for_matrix = "Yes",
      stringsAsFactors = FALSE
    ),
    overwrite = TRUE
  )
  bundle_dir = manifest$Project$bundle_dir[[1]]
  nested_file = file.path(bundle_dir, "28_FeatureMatrixManifest.csv")
  nested = utils::read.csv(nested_file, stringsAsFactors = FALSE,
                           check.names = FALSE)
  nested$FileName[[1]] = "missing_feature_file.csv"
  .categorate_write_csv_file(nested, nested_file)

  validation = validatePlantChemistryAnalysisBundle(bundle_dir)

  expect_equal(validation$Summary$ExportReadyStatus, "fail")
  expect_equal(validation$Summary$ArtifactFailCount, 1)
  expect_equal(validation$Summary$ManifestReferenceFailCount, 1)
  expect_match(validation$ManifestReferences$Message[
    validation$ManifestReferences$Status == "fail"
  ], "missing")
})

test_that("large-run manifests produce retry queues and dry-run plans", {
  manifest = data.frame(
    batch_index = 1:4,
    query_start = c(1, 11, 21, 31),
    query_end = c(10, 20, 30, 40),
    query_count = rep(10, 4),
    query_label = paste("batch", 1:4),
    status = c("completed", "failed", "rate_limited", "not_started"),
    cache_hit_count = c(10, 3, 0, 0),
    request_count = c(0, 7, 4, 0),
    retry_count = c(0, 1, 2, 0),
    error_message = c("", "provider error", "HTTP 503", ""),
    output_file = c("", "", "", ""),
    stringsAsFactors = FALSE
  )

  validation = validatePlantChemistryRunManifest(manifest)
  queue_file = tempfile(fileext = ".csv")
  queue = writePlantChemistryRetryQueue(manifest, queue_file)
  retry = rerunFailedPlantQueries(queue, dry_run = TRUE)

  expect_equal(validation$Summary$batch_count, 4)
  expect_equal(validation$Summary$retry_count, 3)
  expect_equal(nrow(queue), 3)
  expect_true(file.exists(queue_file))
  expect_equal(retry$Plan$status, "planned")
  expect_equal(nrow(retry$RetryRunManifest), 0)
  expect_true(any(grepl("throttle|cooldown|cache|Resume",
                        queue$recommended_action)))
})

test_that("retry runner executes supplied offline batch function", {
  queue = data.frame(
    retry_id = "retry_0001",
    batch_index = 1L,
    query_start = 1L,
    query_end = 2L,
    query_count = 2L,
    query_label = "fixture",
    original_status = "failed",
    retry_reason = "provider_or_batch_error",
    retry_attempt = 1L,
    output_file = NA_character_,
    output_path = NA_character_,
    error_message = "fixture error",
    recommended_action = "rerun",
    created_at = "2026-06-25T00:00:00-0400",
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("retry_run_")
  result = rerunFailedPlantQueries(
    queue,
    dry_run = FALSE,
    out_dir = out_dir,
    overwrite = TRUE,
    runner_fun = function(row) {
      path = file.path(tempdir(), paste0("batch_", row$batch_index, ".rds"))
      saveRDS(list(ok = TRUE), path)
      path
    }
  )

  expect_equal(result$RetryRunManifest$status, "completed")
  expect_true(file.exists(result$RetryRunManifest$output_path))
  expect_true(file.exists(file.path(out_dir, "retry_run_manifest.csv")))
})

test_that("compound identity audit flags review-required chemistry", {
  resolution = data.frame(
    compound_name = c("quercetin", "terpenoid fraction", "beta-pinene"),
    compound_name_clean = c("quercetin", "terpenoid fraction",
                            "beta pinene"),
    query_count = c(5L, 2L, 3L),
    resolved = c(TRUE, FALSE, TRUE),
    CID = c("5280343", NA, "440967"),
    InChIKey = c("REFJWTPEDVJJIY-UHFFFAOYSA-N", NA,
                 "WTARULDDTDQWMU-UHFFFAOYSA-N"),
    SMILES = c("C1=CC(=C(C=C1C2=C(C(=O)C3=C(O2)C=C(C=C3O)O)O)O)O",
               NA, "CC1=CCC(C(C1)C)=C"),
    MolecularFormula = c("C15H10O7", NA, "C10H16"),
    resolution_source = c("pubchem", NA, "pubchem_alias_retry"),
    notes = c("", "mixture label from source", "alias retry"),
    stringsAsFactors = FALSE
  )

  audit = standardizeCompoundIdentityAudit(resolution)
  validation = validateCompoundIdentityAudit(audit)
  template_file = tempfile(fileext = ".csv")
  template = exportCompoundIdentityReviewTemplate(resolution, template_file)

  expect_true(all(c("query_name", "cid", "inchikey_first_block",
                    "review_required", "recommended_action") %in%
                    names(audit)))
  expect_equal(audit$inchikey_first_block[[1]], "REFJWTPEDVJJIY")
  expect_equal(audit$review_required[[2]], "Yes")
  expect_equal(audit$alias_derived_match_flag[[3]], "Yes")
  expect_equal(validation$Summary$identity_count, 3)
  expect_equal(validation$Summary$validation_status, "warn")
  expect_true(file.exists(template_file))
  expect_true(nrow(template) >= 1)
})

test_that("compound identity review wrapper replays explicit decisions", {
  resolution = data.frame(
    compound_name = "unresolved test",
    compound_name_clean = "unresolved test",
    query_count = 1L,
    resolved = FALSE,
    CID = NA_character_,
    InChIKey = NA_character_,
    SMILES = NA_character_,
    MolecularFormula = NA_character_,
    resolution_source = NA_character_,
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
  review = plantCompoundIdentityReviewTable(resolution)
  review$review_decision = "update_identity"
  review$proposed_compound_name = "ethanol"
  review$proposed_cid = "702"
  review$proposed_smiles = "CCO"
  review$proposed_molecular_formula = "C2H6O"
  review$proposed_resolution_source = "manual_identity_review"

  updated = applyCompoundIdentityReview(resolution, review,
                                        reviewer = "unit test")

  expect_true(updated$resolved[[1]])
  expect_equal(updated$CID[[1]], "702")
  expect_equal(updated$SMILES[[1]], "CCO")
  expect_equal(updated$resolution_source[[1]], "manual_identity_review")
})
