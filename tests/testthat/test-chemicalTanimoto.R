tanimoto_fp = function(ch) paste(rep(ch, 154), collapse = "")

tanimoto_request = function(url) {
  decoded = utils::URLdecode(url)
  if (grepl("/pug/compound/inchikey/AAAAAAAAAAAAAA-BBBBBBBBBB-C/cids/JSON$",
            decoded)) {
    return(list(IdentifierList = list(CID = list(101))))
  }
  if (grepl("/pug/compound/inchikey/CCCCCCCCCCCCCC-DDDDDDDDDD-E/cids/JSON$",
            decoded)) {
    return(list(IdentifierList = list(CID = list(102))))
  }
  if (grepl("/pug/compound/inchikey/FFFFFFFFFFFFFF-GGGGGGGGGG-H/cids/JSON$",
            decoded)) {
    return(list(IdentifierList = list(CID = list(103))))
  }
  if (grepl("/property/", decoded)) {
    return(list(PropertyTable = list(Properties = list(
      list(CID = 101,
           Title = "compound alpha",
           MolecularFormula = "C1H1",
           InChIKey = "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
           CanonicalSMILES = "C",
           IsomericSMILES = "C",
           Fingerprint2D = tanimoto_fp("B")),
      list(CID = 102,
           Title = "compound beta",
           MolecularFormula = "C2H2",
           InChIKey = "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
           CanonicalSMILES = "CC",
           IsomericSMILES = "CC",
           Fingerprint2D = tanimoto_fp("C")),
      list(CID = 103,
           Title = "compound gamma",
           MolecularFormula = "C3H3",
           InChIKey = "FFFFFFFFFFFFFF-GGGGGGGGGG-H",
           CanonicalSMILES = "CCC",
           IsomericSMILES = "CCC",
           Fingerprint2D = tanimoto_fp("/"))
    ))))
  }
  list()
}

test_that("chemicalTanimotoSimilarity computes compound and group similarities", {
  input = data.frame(
    species = c("Plant A", "Plant A", "Plant B"),
    compound_name = c("compound alpha", "compound beta", "compound gamma"),
    InChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                 "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                 "FFFFFFFFFFFFFF-GGGGGGGGGG-H"),
    stringsAsFactors = FALSE
  )

  result = chemicalTanimotoSimilarity(
    input,
    group_cols = "species",
    cache = FALSE,
    throttle = 0,
    request_fun = tanimoto_request,
    return_group_compound_pairs = TRUE
  )

  expect_s3_class(result, "uaf_tanimoto_similarity")
  expect_equal(nrow(result$CompoundResolution), 3)
  expect_true(all(result$CompoundResolution$cid_resolution_status ==
                    "inchikey_resolved"))
  expect_equal(nrow(result$PubChemFingerprints), 3)
  expect_equal(nrow(result$CompoundTanimoto), 3)
  expect_true(all(result$CompoundTanimoto$tanimoto >= 0 &
                    result$CompoundTanimoto$tanimoto <= 1))
  expect_equal(nrow(result$GroupPairTanimotoSummary), 1)
  expect_equal(result$GroupPairTanimotoSummary$group_a, "Plant A")
  expect_equal(result$GroupPairTanimotoSummary$group_b, "Plant B")
  expect_equal(result$GroupPairTanimotoSummary$compound_pair_count, 2)
  expect_equal(nrow(result$GroupCompoundTanimoto), 2)
  expect_true(all(c("group_a", "group_b", "compound_id_a",
                    "compound_id_b", "tanimoto") %in%
                    names(result$GroupCompoundTanimoto)))
})

test_that("chemicalTanimotoSimilarity streams large pair tables when out_dir is supplied", {
  input = data.frame(
    species = c("Plant A", "Plant A", "Plant B"),
    compound_name = c("compound alpha", "compound beta", "compound gamma"),
    InChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                 "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                 "FFFFFFFFFFFFFF-GGGGGGGGGG-H"),
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("tanimoto-out-")

  result = chemicalTanimotoSimilarity(
    input,
    group_cols = "species",
    cache = FALSE,
    throttle = 0,
    request_fun = tanimoto_request,
    return_group_compound_pairs = TRUE,
    out_dir = out_dir
  )

  expect_equal(nrow(result$CompoundTanimoto), 0)
  expect_equal(nrow(result$GroupCompoundTanimoto), 0)
  expect_true(file.exists(file.path(out_dir, "compound_pair_tanimoto.csv.gz")))
  expect_true(file.exists(file.path(out_dir,
                                    "group_compound_pair_tanimoto.csv.gz")))
  expect_true(all(c("CompoundTanimoto", "GroupCompoundTanimoto") %in%
                    result$ExportManifest$Table))
})

test_that("plantChemicalTanimotoSimilarity labels plant pair outputs clearly", {
  plant_result = list(
    PlantCompoundOccurrences = data.frame(
      species = c("Plant A", "Plant A", "Plant B"),
      compound_id = c("LTS0001", "LTS0002", "LTS0003"),
      compound_name_clean = c("alpha", "beta", "gamma"),
      compound_name = c("compound alpha", "compound beta", "compound gamma"),
      evidence_tier = "direct_species_database",
      plant_part_group = c("leaf", "root_belowground", "leaf"),
      stringsAsFactors = FALSE
    ),
    CompoundResolution = data.frame(
      compound_name_clean = c("alpha", "beta", "gamma"),
      compound_name = c("compound alpha", "compound beta", "compound gamma"),
      resolved = TRUE,
      CID = NA_character_,
      InChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                   "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                   "FFFFFFFFFFFFFF-GGGGGGGGGG-H"),
      SMILES = c("C", "CC", "CCC"),
      MolecularFormula = c("C1H1", "C2H2", "C3H3"),
      resolution_source = "test_fixture",
      stringsAsFactors = FALSE
    )
  )

  result = plantChemicalTanimotoSimilarity(
    plant_result,
    cache = FALSE,
    throttle = 0,
    request_fun = tanimoto_request,
    return_group_compound_pairs = TRUE
  )

  expect_s3_class(result, "uaf_plant_tanimoto_similarity")
  expect_true("species_a" %in% names(result$PlantPairTanimotoSummary))
  expect_true("species_b" %in% names(result$PlantPairTanimotoSummary))
  expect_equal(result$PlantPairTanimotoSummary$species_a, "Plant A")
  expect_equal(result$PlantPairTanimotoSummary$species_b, "Plant B")
  expect_true("species_compound_pair_id" %in%
                names(result$PlantCompoundTanimoto))
  expect_equal(nrow(result$PlantCompoundMembership), 3)
  expect_true("source_compound_id" %in% names(result$PlantCompoundMembership))
  expect_true(all(grepl("^[A-Z]{14}_", result$PlantCompoundMembership$compound_id)))
})
