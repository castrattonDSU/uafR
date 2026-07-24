categorate_batch_fixture = function(query, cid, property_cids = cid,
                                    validation_status = "pass") {
  data.frame_safe = function(...) data.frame(..., stringsAsFactors = FALSE)
  list(
    PubChemIdentity = data.frame_safe(
      Query = query,
      CID = cid,
      MatchStatus = "resolved_cid"
    ),
    PubChemProperties = data.frame_safe(
      Query = query[match(property_cids, cid)],
      CID = property_cids,
      MolecularFormula = paste0("C", seq_along(property_cids), "H2")
    ),
    ChemicalTraits = data.frame_safe(
      Query = query,
      CID = cid,
      TraitType = "property",
      TraitGroup = "formula_present",
      TraitValue = "yes",
      MatrixKey = "property__formula_present__yes",
      Confidence = "high"
    ),
    ChemicalTraitMatrix = data.frame_safe(
      Query = query,
      CID = cid,
      property__formula_present__yes = 1
    ),
    ChemicalTraitEvidence = data.frame_safe(
      Query = query,
      CID = cid,
      EvidenceType = "trait",
      AnalysisKey = "property__formula_present__yes",
      EvidenceText = "Fixture evidence"
    ),
    DerivedGroups = data.frame_safe(
      Query = query,
      CID = cid,
      natural_product_classes = c(
        "phenolics, quoted \"class\" | newline normalized",
        "terpenoids"
      )[seq_along(query)],
      metabolic_context = "specialized metabolite",
      ecological_context = "fixture context",
      lipophilicity_bin = "not_reported",
      volatility_proxy = "not_reported",
      oxygenated = TRUE,
      nitrogenous = FALSE,
      sulfur_containing = FALSE
    ),
    SourceCoverage = data.frame_safe(
      Query = query,
      Source = "PubChemProperties",
      RecordCount = as.integer(cid %in% property_cids)
    ),
    ValidationSummary = data.frame_safe(
      Status = validation_status,
      IssueCount = ifelse(validation_status == "pass", 0L, 1L)
    ),
    ValidationIssues = if (validation_status == "pass") {
      data.frame()
    } else {
      data.frame_safe(
        Table = "Fixture",
        Severity = "warning",
        Issue = "Fixture warning"
      )
    },
    DataDictionary = data.frame_safe(
      Table = "ChemicalTraits",
      Column = "Query",
      Type = "character"
    ),
    Provenance = data.frame_safe(
      Function = "categorate_batch_fixture",
      RetrievedAt = "2026-06-23T00:00:00-0400"
    )
  )
}

write_categorate_batch_fixtures = function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    categorate_batch_fixture(
      query = c("cid:1", "cid:2"),
      cid = c(1L, 2L)
    ),
    file.path(path, "categorate_batch_0001.rds")
  )
  saveRDS(
    categorate_batch_fixture(
      query = c("cid:3", "cid:4"),
      cid = c(3L, 4L),
      property_cids = 3L,
      validation_status = "warn"
    ),
    file.path(path, "categorate_batch_0002.rds")
  )
  saveRDS(
    simpleError("fixture batch failed"),
    file.path(path, "categorate_batch_0003.rds")
  )
  invisible(path)
}

test_that("categorate batch reader and summary flag quality issues", {
  batch_dir = write_categorate_batch_fixtures(tempfile("categorate_batches_"))

  batches = readCategorateBatchDirectory(batch_dir, min_property_ratio = 0.9)
  summary = batches$BatchSummary

  expect_s3_class(batches, "uaf_categorate_batch_directory")
  expect_equal(nrow(summary), 3)
  expect_equal(summary$Status, c("ok", "incomplete", "error"))
  expect_equal(summary$ResolvedCIDCount[1:2], c(2L, 2L))
  expect_equal(summary$PubChemPropertyRows[1:2], c(2L, 1L))
  expect_equal(summary$PubChemPropertyRatio[2], 0.5)
  expect_match(summary$QualityIssue[2], "Only 1 PubChem property rows")
  expect_match(summary$Error[3], "fixture batch failed")
})

test_that("categorate batch combiner binds successful batches with provenance", {
  batch_dir = write_categorate_batch_fixtures(tempfile("categorate_batches_"))

  combined = combineCategorateTables(
    batch_dir,
    tables = c("ChemicalTraits", "PubChemProperties")
  )

  expect_equal(names(combined), c("ChemicalTraits", "PubChemProperties"))
  expect_equal(nrow(combined$ChemicalTraits), 4)
  expect_equal(nrow(combined$PubChemProperties), 3)
  expect_true(all(c("CategorateBatchIndex", "CategorateBatchFile",
                    "Query", "CID") %in% names(combined$ChemicalTraits)))
  expect_false(any(combined$ChemicalTraits$CategorateBatchIndex == 3))
  expect_equal(sort(unique(combined$ChemicalTraits$CategorateBatchIndex)),
               c(1L, 2L))
})

test_that("external file references are portable and checksum-aware", {
  artifact = tempfile("uafr_external_artifact_", fileext = ".csv")
  writeLines(c("a,b", "1,2"), artifact, useBytes = TRUE)

  refs = .categorate_analysis_file_references(
    c(pairwise_tanimoto = artifact,
      missing_artifact = file.path(tempdir(), "missing_pairwise.csv.gz"))
  )

  expect_equal(refs$Reference,
               c("pairwise_tanimoto", "missing_artifact"))
  expect_equal(refs$Path,
               c(basename(artifact), "missing_pairwise.csv.gz"))
  expect_true(all(refs$PathType == "external_reference_basename"))
  expect_equal(refs$OriginalPathWasAbsolute, c("Yes", "Yes"))
  expect_equal(refs$ReferenceStatus,
               c("available_at_export", "missing_at_export"))
  expect_equal(refs$ChecksumStatus, c("computed", "missing_at_export"))
  expect_equal(refs$ArtifactChecksum[[1]], unname(tools::md5sum(artifact)))
  expect_true(is.na(refs$ArtifactChecksum[[2]]))
  expect_false(any(grepl("^/", refs$Path)))
})

test_that("plant chemistry analysis bundle exports combined batch tables", {
  batch_dir = write_categorate_batch_fixtures(tempfile("categorate_batches_"))
  membership = data.frame(
    species = c("Plant alpha", "Plant beta"),
    compound_id = c("cid_1", "cid_3"),
    compound_name = c("compound one", "compound three"),
    evidence_tier = "direct_species_database",
    stringsAsFactors = FALSE
  )
  species_pairs = data.frame(
    species_a = "Plant alpha",
    species_b = "Plant beta",
    mean_tanimoto = 0.75,
    compound_pair_count = 1L,
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("plant_chemistry_bundle_")

  manifest = exportPlantChemistryAnalysisBundle(
    categorate_batches = batch_dir,
    path = out_dir,
    plant_membership = membership,
    species_pair_tanimoto = species_pairs,
    file_references = c(
      plant_compound_pairs = "large_plant_compound_pair_tanimoto.csv.gz"
    ),
    tables = c("ChemicalTraits", "PubChemProperties"),
    format = "csv",
    overwrite = TRUE
  )

  expect_true(dir.exists(out_dir))
  expect_true("ExportManifest" %in% manifest$Table)
  expect_true("DataDictionary" %in% manifest$Table)
  expect_true(all(c("BatchSummary", "PlantCompoundMembership",
                    "PlantPairTanimotoSummary", "FileReferences",
                    "ChemicalTraits", "PubChemProperties") %in%
                    manifest$Table))
  traits_file = file.path(
    out_dir,
    manifest$FileName[manifest$Table == "ChemicalTraits"]
  )
  batch_file = file.path(
    out_dir,
    manifest$FileName[manifest$Table == "BatchSummary"]
  )
  refs_file = file.path(
    out_dir,
    manifest$FileName[manifest$Table == "FileReferences"]
  )
  exported_traits = utils::read.csv(traits_file, stringsAsFactors = FALSE,
                                    check.names = FALSE)
  exported_summary = utils::read.csv(batch_file, stringsAsFactors = FALSE,
                                     check.names = FALSE)
  exported_refs = utils::read.csv(refs_file, stringsAsFactors = FALSE,
                                  check.names = FALSE)

  expect_equal(nrow(exported_traits), 4)
  expect_equal(nrow(exported_summary), 3)
  expect_equal(exported_summary$Status, c("ok", "incomplete", "error"))
  expect_equal(exported_refs$Reference, "plant_compound_pairs")
  expect_false(exported_refs$Exists)
  expect_equal(exported_refs$Path,
               "large_plant_compound_pair_tanimoto.csv.gz")
  expect_equal(exported_refs$PathType, "external_reference_basename")
  expect_equal(exported_refs$ChecksumStatus, "missing_at_export")
  expect_equal(exported_refs$ReferenceStatus, "missing_at_export")

  enriched_file = file.path(out_dir, "03b_PlantCompoundMembershipEnriched.csv")
  species_file = file.path(out_dir, "14_PlantChemistrySummary.csv")
  missing_file = file.path(out_dir, "15_PlantChemistryMissingSpecies.csv")
  validation_file = file.path(out_dir, "12b_ValidationOverview.csv")

  expect_true(file.exists(enriched_file))
  expect_true(file.exists(species_file))
  expect_true(file.exists(missing_file))
  expect_true(file.exists(validation_file))

  enriched = utils::read.csv(enriched_file, stringsAsFactors = FALSE,
                             check.names = FALSE)
  expect_true(all(c("taxonomy_family_status", "has_fingerprint",
                    "comparison_scope", "comparison_group",
                    "comparable_for_matrix") %in% names(enriched)))
  validation = validatePlantChemistryAnalysisBundle(out_dir)
  expect_equal(validation$Summary$ExportReadyStatus, "pass")
  expect_equal(validation$Summary$ArtifactFailCount, 0)
})

test_that("streamed bundle CSV writer preserves parseable quoted fields", {
  batch_dir = write_categorate_batch_fixtures(tempfile("categorate_batches_"))
  out_dir = tempfile("plant_chemistry_bundle_csv_")

  manifest = exportPlantChemistryAnalysisBundle(
    categorate_batches = batch_dir,
    path = out_dir,
    plant_membership = data.frame(
      species = "Plant alpha",
      compound_id = "cid_1",
      compound_name = "compound one",
      source_database = "fixture",
      evidence_tier = "direct_species_database",
      stringsAsFactors = FALSE
    ),
    tables = c("DerivedGroups"),
    format = "csv",
    overwrite = TRUE
  )

  derived_file = file.path(
    out_dir,
    manifest$FileName[manifest$Table == "DerivedGroups"]
  )
  fields = utils::count.fields(derived_file, sep = ",", quote = "\"",
                               blank.lines.skip = FALSE)
  expect_true(length(fields) > 1)
  expect_true(all(fields == fields[[1]]))

  derived = utils::read.csv(derived_file, stringsAsFactors = FALSE,
                            check.names = FALSE)
  expect_equal(nrow(derived), 4)
  expect_true("natural_product_classes" %in% names(derived))
  expect_equal(
    validatePlantChemistryAnalysisBundle(out_dir)$Summary$ExportReadyStatus,
    "pass"
  )
})
