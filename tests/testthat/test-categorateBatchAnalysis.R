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
  expect_equal(manifest$Table[[1]], "ExportManifest")
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
})
