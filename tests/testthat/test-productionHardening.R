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
  plan = planPlantChemistryRun(
    plants = paste("Plant", seq_len(100)),
    compounds = data.frame(compound_id = paste0("c", seq_len(25))),
    sources = c("lotus", "pubmed"),
    cache_dir = cache_dir,
    lotus_index = "local_lotus_lookup",
    write_full_pairwise = TRUE
  )

  expect_equal(nrow(inventory), 2)
  expect_true("PubChem" %in% summary$provider)
  expect_equal(plan$Summary$species_count, 100)
  expect_equal(plan$OutputEstimates$estimated_unique_compounds, 25)
  expect_true(any(plan$ProviderPlan$recommended_mode == "local_index"))
})

test_that("chemistry classification overrides are reproducible", {
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
  expect_equal(out$comparison_group[[1]], "terpenoids")
  expect_equal(out$classification_source[[1]], "user_override")
  expect_equal(out$comparable_for_matrix[[1]], "Yes")
  expect_equal(out$comparison_scope[[2]], "unknown")
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
    modes = c("count", "binary", "fraction", "confidence")
  )

  expect_equal(nrow(direct), 2)
  expect_equal(nrow(comparable), 2)
  expect_equal(nrow(review), 1)
  expect_true(all(c("EvidenceGradeCountMatrix",
                    "PlantPartCountMatrix",
                    "ComparisonGroupConfidenceMatrix") %in% names(features)))
  expect_true("confidence" %in% features$Manifest$Mode)
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
                    "workflow_name", "workflow_parameters") %in%
                    names(bundle_manifest)))
  expect_true("FeatureMatrixManifest" %in% bundle_manifest$Table)
  expect_equal(
    validatePlantChemistryAnalysisBundle(
      manifest$Project$bundle_dir[[1]]
    )$Summary$ExportReadyStatus,
    "pass"
  )
})
