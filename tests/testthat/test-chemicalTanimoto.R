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

test_that("Tanimoto preparation collapses verified aliases and rejects identity mismatches", {
  alpha_key = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
  beta_key = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"
  wrong_key = "FFFFFFFFFFFFFF-GGGGGGGGGG-H"
  request_fun = function(url) {
    decoded = utils::URLdecode(url)
    if (grepl(paste0("/inchikey/", alpha_key, "/cids/JSON$"), decoded) ||
        grepl(paste0("/inchikey/", wrong_key, "/cids/JSON$"), decoded)) {
      return(list(IdentifierList = list(CID = list(101))))
    }
    if (grepl(paste0("/inchikey/", beta_key, "/cids/JSON$"), decoded)) {
      return(list(IdentifierList = list(CID = list(102))))
    }
    if (grepl("/property/", decoded)) {
      return(list(PropertyTable = list(Properties = list(
        list(CID = 101, Title = "alpha", MolecularFormula = "C1H1",
             InChIKey = alpha_key, CanonicalSMILES = "C",
             IsomericSMILES = "C", Fingerprint2D = tanimoto_fp("B")),
        list(CID = 102, Title = "beta", MolecularFormula = "C2H2",
             InChIKey = beta_key, CanonicalSMILES = "CC",
             IsomericSMILES = "CC", Fingerprint2D = tanimoto_fp("C"))
      ))))
    }
    list()
  }
  input = data.frame(
    species = c("Plant A", "Plant A", "Plant B", "Plant B"),
    source_id = c("alpha_primary", "alpha_alias", "beta_primary",
                  "bad_cached_identity"),
    compound_name = c("alpha", "alpha synonym", "beta", "wrong identity"),
    InChIKey = c(alpha_key, alpha_key, beta_key, wrong_key),
    source_record_id = c("r1", "r2", "r3", "r4"),
    source_database = "fixture",
    evidence_tier = "direct_species_database",
    stringsAsFactors = FALSE
  )

  result = chemicalTanimotoSimilarity(
    input,
    compound_id_col = "source_id",
    group_cols = "species",
    cache = FALSE,
    throttle = 0,
    request_fun = request_fun
  )

  expect_equal(nrow(result$CompoundInput), 4L)
  expect_equal(nrow(result$PubChemFingerprints), 4L)
  expect_equal(nrow(result$CompoundIdentityMap), 4L)
  expect_equal(nrow(result$CanonicalCompounds), 2L)
  expect_equal(nrow(result$ExcludedCompoundIdentities), 1L)
  expect_equal(result$ExcludedCompoundIdentities$input_compound_id,
               "bad_cached_identity")
  expect_equal(result$ExcludedCompoundIdentities$identity_match_status,
               "input_pubchem_inchikey_mismatch")
  expect_equal(nrow(result$GroupCompoundMembershipEvidence), 3L)
  expect_equal(nrow(result$GroupCompoundMembership), 2L)
  alpha = result$GroupCompoundMembership[
    result$GroupCompoundMembership$group_id == "Plant A", , drop = FALSE
  ]
  expect_equal(alpha$collapsed_membership_row_count, 2L)
  expect_match(alpha$tanimoto_input_compound_ids, "alpha_alias")
  expect_match(alpha$tanimoto_input_compound_ids, "alpha_primary")
  expect_equal(nrow(result$CompoundTanimoto), 1L)
  expect_equal(result$GroupPairTanimotoSummary$compound_pair_count, 1L)
})

test_that("preparePlantTanimotoInput uses exact source records and writes a safe handoff", {
  alpha_key = "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
  quercetin_a = "CCCCCCCCCCCCCC-DDDDDDDDDD-E"
  quercetin_b = "EEEEEEEEEEEEEE-FFFFFFFFFF-G"
  review_key = "GGGGGGGGGGGGGG-HHHHHHHHHH-I"
  low_key = "IIIIIIIIIIIIII-JJJJJJJJJJ-K"
  not_ready_key = "KKKKKKKKKKKKKK-LLLLLLLLLL-M"
  occurrences = data.frame(
    species = c("Plant A", "Plant A", "Plant B", "Plant C", "Plant D",
                "Plant E", "Plant F", "Plant G"),
    compound_id = c("LTS_ALPHA", "LTS_ALPHA", "LTS_QA", "LTS_QB",
                    "LTS_CLASS", "LTS_MISSING", "LTS_LOW", "LTS_NOT_READY"),
    compound_name = c("alpha", "alpha synonym", "quercetin", "quercetin",
                      "flavonoids", "unknown compound", "low evidence",
                      "not ready"),
    compound_name_clean = c("alpha", "alpha", "quercetin", "quercetin",
                            "flavonoids", "unknown_compound", "low_evidence",
                            "not_ready"),
    source_database = "LOTUS",
    source_record_id = c("LTS_ALPHA", "LTS_ALPHA", "LTS_QA", "LTS_QB",
                         "LTS_CLASS", "LTS_MISSING", "LTS_LOW",
                         "LTS_NOT_READY"),
    occurrence_status = "direct_reported",
    analysis_ready = c(rep("Yes", 7), NA_character_),
    confidence = c(rep("high", 6), "low", "high"),
    CID = "67761025",
    Fingerprint2D = "legacy_fingerprint_must_not_be_reused",
    stringsAsFactors = FALSE
  )
  source_identity = data.frame(
    compound_name = c("alpha", "quercetin", "quercetin", "flavonoids",
                      "low evidence", "not ready"),
    compound_name_clean = c("alpha", "quercetin", "quercetin", "flavonoids",
                            "low_evidence", "not_ready"),
    source_database = "LOTUS",
    source_record_id = c("LTS_ALPHA", "LTS_QA", "LTS_QB", "LTS_CLASS",
                         "LTS_LOW", "LTS_NOT_READY"),
    source_compound_id = c("LTS_ALPHA", "LTS_QA", "LTS_QB", "LTS_CLASS",
                           "LTS_LOW", "LTS_NOT_READY"),
    source_compound_id_type = "LOTUS",
    CID = NA_character_,
    InChIKey = c(alpha_key, quercetin_a, quercetin_b, review_key, low_key,
                  not_ready_key),
    SMILES = c("C", "CC", "CCC", "CCCC", "CCCCC", "CCCCCC"),
    MolecularFormula = c("CH4", "C2H6", "C3H8", "C4H10", "C5H12",
                         "C6H14"),
    evidence_url = "https://lotus.test/download",
    evidence_text = "source fixture",
    identity_status = c("source_structure_unique",
                        "source_structure_ambiguous",
                        "source_structure_ambiguous",
                        rep("source_structure_unique", 3)),
    identity_note = "",
    stringsAsFactors = FALSE
  )
  review = data.frame(
    compound_name_clean = c("quercetin", "flavonoids"),
    identity_issue_type = c("resolved_source_structure_ambiguous",
                            "resolved_broad_class_label"),
    review_decision = NA_character_,
    stringsAsFactors = FALSE
  )
  input = list(
    PlantCompoundOccurrences = occurrences,
    CompoundResolution = data.frame(),
    SourceCompoundIdentity = source_identity,
    CompoundIdentityReview = review
  )
  out_dir = tempfile("tanimoto-handoff-")

  prepared = preparePlantTanimotoInput(
    input,
    out_dir = out_dir,
    strict = TRUE
  )

  expect_s3_class(prepared, "uaf_plant_tanimoto_input")
  expect_equal(nrow(prepared$PlantCompoundMembership), 3L)
  expect_equal(nrow(prepared$CompoundInput), 3L)
  expect_equal(nrow(prepared$DuplicateAudit), 1L)
  expect_equal(prepared$DuplicateAudit$evidence_row_count, 2L)
  expect_equal(nrow(prepared$NameStructureAudit), 1L)
  expect_identical(prepared$NameStructureAudit$audit_scope, "global")
  expect_identical(prepared$NameStructureAudit$compound_name_clean,
                   "quercetin")
  expect_equal(prepared$NameStructureAudit$structure_count, 2L)
  expect_equal(prepared$NameStructureAudit$group_count, 2L)
  expect_match(prepared$NameStructureAudit$groups, "Plant B")
  expect_match(prepared$NameStructureAudit$groups, "Plant C")
  expect_true(all(prepared$PlantCompoundMembership$source_identity_exact))
  expect_true(all(is.na(prepared$PlantCompoundMembership$CID)))
  expect_true(all(grepl("67761025",
                        prepared$PlantCompoundMembership$reported_CIDs,
                        fixed = TRUE)))
  expect_false(any(grepl("fingerprint", names(prepared$EvidenceRows),
                         ignore.case = TRUE)))
  expect_equal(prepared$Summary$ignored_input_fingerprint_column_count, 1L)
  expect_true(any(grepl("identity_review_required",
                        prepared$ExcludedRows$preparation_exclusion_reason,
                        fixed = TRUE)))
  expect_true(any(grepl("below_minimum_confidence",
                        prepared$ExcludedRows$preparation_exclusion_reason,
                        fixed = TRUE)))
  expect_true(any(grepl("not_analysis_ready",
                        prepared$ExcludedRows$preparation_exclusion_reason,
                        fixed = TRUE)))
  expect_true(any(grepl("missing_valid_inchikey_or_cid",
                        prepared$ExcludedRows$preparation_exclusion_reason,
                        fixed = TRUE)))
  expect_true(all(prepared$ValidationSummary$status == "pass"))

  expect_equal(nrow(prepared$ExportManifest), 9L)
  for (i in seq_len(nrow(prepared$ExportManifest))) {
    path = file.path(out_dir, prepared$ExportManifest$file[[i]])
    expect_true(file.exists(path))
    expect_identical(unname(tools::md5sum(path)[[1]]),
                     prepared$ExportManifest$md5[[i]])
  }
  expect_true(file.exists(file.path(out_dir,
                                    "tanimoto_input_export_manifest.csv")))
  manifest = jsonlite::read_json(
    file.path(out_dir, "tanimoto_input_manifest.json"),
    simplifyVector = TRUE
  )
  expect_identical(manifest$schema_version, "1.0.0")
  expect_false(any(grepl("pair_tanimoto|compound_pair",
                         list.files(out_dir), ignore.case = TRUE)))
})

test_that("preparePlantTanimotoInput does not write a strict invalid handoff", {
  out_dir = tempfile("invalid-tanimoto-handoff-")
  input = data.frame(
    species = "Plant A",
    compound_name = "compound alpha",
    InChIKey = "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
    occurrence_status = "direct_reported",
    analysis_ready = "Yes",
    confidence = "high",
    stringsAsFactors = FALSE
  )

  expect_error(
    preparePlantTanimotoInput(input, out_dir = out_dir, strict = TRUE),
    "at_least_two_structures"
  )
  expect_false(dir.exists(out_dir))
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

test_that("Tanimoto smoke selection is deterministic and plant-balanced", {
  input = data.frame(
    species = c("Plant C", "Plant A", "Plant B", "Plant A", "Plant B"),
    compound_id = c("c5", "c1", "c3", "c2", "c4"),
    compound_name = paste("compound", 1:5),
    stringsAsFactors = FALSE
  )
  selected_a = .tanimoto_smoke_membership(input, n = 3)
  selected_b = .tanimoto_smoke_membership(input[5:1, ], n = 3)

  expect_identical(attr(selected_a, "selected_compound_ids"),
                   c("c1", "c3", "c5"))
  expect_identical(attr(selected_a, "selected_compound_ids"),
                   attr(selected_b, "selected_compound_ids"))
  expect_setequal(selected_a$species, c("Plant A", "Plant B", "Plant C"))
})

test_that("chemicalTanimotoSimilarity shards pair outputs and reports progress", {
  input = data.frame(
    species = c("Plant A", "Plant A", "Plant B"),
    compound_name = c("compound alpha", "compound beta", "compound gamma"),
    InChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                 "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
                 "FFFFFFFFFFFFFF-GGGGGGGGGG-H"),
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("tanimoto-shards-")
  events = list()

  result = chemicalTanimotoSimilarity(
    input,
    group_cols = "species",
    cache = FALSE,
    throttle = 0,
    request_fun = tanimoto_request,
    out_dir = out_dir,
    pair_shard_rows = 2,
    pair_block_size = 1,
    progress_every = 1,
    progress_fun = function(event) events[[length(events) + 1L]] <<- event
  )

  pair_manifest = result$ExportManifest[
    result$ExportManifest$Table == "CompoundTanimoto", , drop = FALSE
  ]
  expect_equal(nrow(pair_manifest), 2L)
  expect_equal(sum(as.numeric(pair_manifest$Rows)), 3)
  expect_true(all(file.exists(pair_manifest$Path)))
  expect_false(any(grepl("[.]partial$", list.files(out_dir))))
  event_names = vapply(events, function(event) event$event, character(1))
  expect_true("identities_completed" %in% event_names)
  expect_true("property_chunk_completed" %in% event_names)
  expect_true("pair_blocks_completed" %in% event_names)
})

test_that("incomplete pair writers never expose a completed gzip path", {
  path = tempfile("incomplete-pairs-", fileext = ".csv.gz")
  writer = .tanimoto_pair_writer(
    path = path,
    columns = c("a", "b"),
    expected_rows = 2,
    shard_rows = Inf
  )
  writer$write(data.frame(a = 1, b = 2))

  expect_error(writer$finalize(), "produced 1 rows; expected 2")
  writer$abort()
  expect_false(file.exists(path))
  expect_true(file.exists(paste0(path, ".partial")))
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

test_that("Tanimoto preparation preserves conservative chemistry classes", {
  rows = data.frame(
    species = c("Plant A", "Plant B"),
    compound_name = c("alpha", "beta"),
    compound_name_clean = c("alpha", "beta"),
    InChIKey = c("AAAAAAAAAAAAAA-BBBBBBBBBB-C",
                  "CCCCCCCCCCCCCC-DDDDDDDDDD-E"),
    CID = NA_character_,
    occurrence_status = "direct_reported",
    analysis_ready = "Yes",
    confidence = "high",
    comparison_scope = "specialized_metabolites",
    comparison_group = "terpenoids",
    comparison_subgroup = "monoterpenoids",
    metabolism_domain = "specialized_metabolism",
    biosynthetic_family = "terpenoid",
    chemical_behavior = "volatile",
    comparability_confidence = "high",
    comparability_basis = "source-backed fixture",
    classification_source = "PubChemClassifications",
    comparable_for_matrix = "Yes",
    comparison_caveat =
      "Source-backed classification; sample occurrence is not confirmed.",
    stringsAsFactors = FALSE
  )

  prepared = preparePlantTanimotoInput(rows, strict = TRUE)

  expect_equal(prepared$PlantCompoundMembership$comparison_scope,
               rep("specialized_metabolites", 2))
  expect_equal(prepared$PlantCompoundMembership$comparison_group,
               rep("terpenoids", 2))
  expect_equal(prepared$PlantCompoundMembership$comparable_for_matrix,
               rep("Yes", 2))
  expect_true(all(c("comparability_confidence", "classification_sources",
                    "comparison_caveats") %in%
                    names(prepared$PlantCompoundMembership)))
})

test_that("decoded fingerprints produce exact comparable summaries", {
  bits = rbind(
    c1 = c(TRUE, TRUE, FALSE, FALSE),
    c2 = c(TRUE, FALSE, TRUE, FALSE),
    c3 = c(FALSE, FALSE, TRUE, TRUE)
  )
  membership = data.frame(
    group_id = c("Plant A", "Plant B", "Plant B"),
    compound_id = c("c1", "c2", "c3"),
    comparison_scope = c("specialized_metabolites",
                         "specialized_metabolites", "primary_metabolites"),
    comparison_group = c("terpenoids", "terpenoids", "carbohydrates"),
    comparable_for_matrix = "Yes",
    stringsAsFactors = FALSE
  )
  compounds = data.frame(
    compound_id = c("c1", "c2", "c3"),
    compound_name = c("alpha", "beta", "gamma"),
    stringsAsFactors = FALSE
  )

  scope = .tanimoto_comparable_pair_summaries(
    bits, membership, compounds, "species", c(0.5, 0.7), 5,
    field = "comparison_scope", type = "scope"
  )
  group = .tanimoto_comparable_pair_summaries(
    bits, membership, compounds, "species", c(0.5, 0.7), 5,
    field = "comparison_group", type = "group"
  )

  expect_equal(nrow(scope), 1L)
  expect_equal(scope$comparison_scope, "specialized_metabolites")
  expect_equal(scope$compound_pair_count, 1L)
  expect_equal(scope$mean_tanimoto, 1 / 3, tolerance = 1e-6)
  expect_equal(group$comparison_group, "terpenoids")
})
