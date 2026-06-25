test_that("plant metadata normalization preserves supplied taxonomy status", {
  meta = data.frame(
    Scientific.Name = c("Salix nigra", "Zea mays"),
    Accepted.Name = c("Salix nigra", ""),
    Plant.Family = c("Salicaceae", "Poaceae"),
    Project.Role = c("candidate", "reference"),
    stringsAsFactors = FALSE
  )

  out = standardizePlantMetadata(meta)

  expect_equal(out$species, c("Salix nigra", "Zea mays"))
  expect_equal(out$genus, c("Salix", "Zea"))
  expect_equal(out$family_status, c("family_supplied", "family_supplied"))
  expect_equal(out$accepted_name_status,
               c("accepted_name_supplied", "accepted_name_not_supplied"))
  expect_equal(out$role, c("candidate", "reference"))
})

test_that("plant occurrence evidence grades stay conservative", {
  occ = data.frame(
    species = c("A species", "B species", "C species", "D species"),
    compound_name = c("alpha", "beta", "gamma", ""),
    compound_name_clean = c("alpha", "beta", "gamma", ""),
    compound_id = c("cid_1", "cid_2", "cid_3", ""),
    source_database = c("LOTUS", "PubMed", "PubTator", "manual"),
    source_record_id = c("L1", "", "", ""),
    evidence_tier = c("direct_species_database",
                      "genus_database_fallback",
                      "direct_species_pubtator_candidate",
                      "unresolved"),
    matched_rank = c("species", "genus", "species", "species"),
    SMILES = c("CCO", "CCC", NA, NA),
    comparable_for_matrix = c("Yes", "Yes", "No", "No"),
    comparison_scope = c("specialized_metabolites",
                         "specialized_metabolites", "unknown", "unknown"),
    comparison_group = c("phenolics", "terpenoids", "unknown", "unknown"),
    stringsAsFactors = FALSE
  )
  review = data.frame(row_id = 2L, curation_decision = "exclude",
                      stringsAsFactors = FALSE)

  out = plantOccurrenceEvidenceGrade(occ, review_table = review)

  expect_equal(out$evidence_grade[[1]], "direct_species_database_record")
  expect_equal(out$evidence_grade[[2]], "excluded_by_review")
  expect_equal(out$evidence_grade[[3]], "pubtator_pubmed_candidate_only")
  expect_equal(out$evidence_grade[[4]], "unresolved_or_review_required")
  expect_equal(out$source_backed[[1]], "Yes")
  expect_equal(out$review_required[[3]], "Yes")
})

test_that("comparable Tanimoto summaries filter by matching scope and group", {
  membership = data.frame(
    species = c("Plant A", "Plant B", "Plant A", "Plant B", "Plant B"),
    compound_id = c("c1", "c2", "c3", "c4", "c5"),
    compound_name = paste0("compound", 1:5),
    comparison_scope = c("specialized_metabolites",
                         "specialized_metabolites",
                         "primary_metabolites", "primary_metabolites",
                         "unknown"),
    comparison_group = c("terpenoids", "terpenoids",
                         "amino_acids", "amino_acids", "unknown"),
    comparable_for_matrix = c("Yes", "Yes", "Yes", "Yes", "No"),
    source_database = "fixture",
    evidence_tier = "direct_species_database",
    stringsAsFactors = FALSE
  )
  pairs = data.frame(
    species_a = "Plant A",
    compound_id_a = c("c1", "c1", "c3", "c3"),
    compound_name_a = c("compound1", "compound1", "compound3", "compound3"),
    species_b = "Plant B",
    compound_id_b = c("c2", "c4", "c4", "c5"),
    compound_name_b = c("compound2", "compound4", "compound4", "compound5"),
    tanimoto = c(0.8, 0.2, 0.6, 0.9),
    stringsAsFactors = FALSE
  )

  out = plantComparableTanimotoSummary(
    plant_compound_pair_tanimoto = pairs,
    membership = membership,
    thresholds = c(0.25, 0.75)
  )

  expect_equal(sort(out$ScopeFiltered$comparison_scope),
               c("primary_metabolites", "specialized_metabolites"))
  expect_equal(sort(out$GroupFiltered$comparison_group),
               c("amino_acids", "terpenoids"))
  expect_false(any(out$ScopeFiltered$comparison_value == "unknown"))
  expect_equal(out$ScopeFiltered$compound_pair_count, c(1L, 1L))
  expect_true(all(c("compound_pair_count_ge_0_25",
                    "compound_pair_count_ge_0_75") %in%
                    names(out$ScopeFiltered)))
})

test_that("feature exports produce stable species matrices and metadata", {
  membership = data.frame(
    species = c("Plant A", "Plant A", "Plant B"),
    compound_id = c("c1", "c2", "c3"),
    compound_name = c("alpha", "beta", "gamma"),
    comparison_scope = c("specialized_metabolites",
                         "specialized_metabolites", "primary_metabolites"),
    comparison_group = c("terpenoids", "phenolics", "amino_acids"),
    comparable_for_matrix = "Yes",
    source_database = c("LOTUS", "LOTUS", "PubMed"),
    evidence_tier = "direct_species_database",
    SMILES = "CCO",
    stringsAsFactors = FALSE
  )

  out = exportPlantChemistryFeatureSet(membership)

  expect_true(all(c("species_id", "species") %in%
                    names(out$ComparisonGroupCountMatrix)))
  expect_equal(out$SpeciesMetadata$compound_count,
               c(2L, 1L))
  expect_equal(out$SpeciesMetadata$species_id,
               c("species__plant_a", "species__plant_b"))
})

test_that("offline plant chemistry project runner finalizes a reusable bundle", {
  out_dir = tempfile("uafr_project_runner_")
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
  pairs = data.frame(
    species_a = "Plant A",
    compound_id_a = "c1",
    species_b = "Plant B",
    compound_id_b = "c2",
    tanimoto = 0.75,
    stringsAsFactors = FALSE
  )

  manifest = runPlantChemistryProject(
    plant_list = c("Plant A", "Plant B", "Plant C"),
    output_dir = out_dir,
    plant_compounds = plant_compounds,
    plant_compound_pair_tanimoto = pairs,
    project_id = "test_project",
    overwrite = TRUE
  )

  bundle_dir = manifest$Project$bundle_dir[[1]]
  expect_true(file.exists(file.path(bundle_dir, "16_EvidenceGradeSummary.csv")))
  expect_true(file.exists(file.path(bundle_dir,
                                    "18_FeatureComparisonGroupCountMatrix.csv")))
  expect_true(file.exists(file.path(bundle_dir,
                                    "22_ComparableScopeTanimotoSummary.csv")))
  validation = validatePlantChemistryAnalysisBundle(bundle_dir)
  expect_equal(validation$Summary$ExportReadyStatus, "pass")
  missing = utils::read.csv(file.path(bundle_dir,
                                      "15_PlantChemistryMissingSpecies.csv"),
                            stringsAsFactors = FALSE, check.names = FALSE)
  expect_equal(missing$species, "Plant C")
})
