mock_ainsect_pubchem_profile = function(compounds,
                                        profile = "ms",
                                        cache = TRUE,
                                        cache_dir = NULL,
                                        throttle = 0.2,
                                        ...) {
  stopifnot(profile == "ms")
  identity = data.frame(
    Query = compounds,
    CID = c(Limonene = 22311, Linalool = 6549, `Unknown EO` = NA_integer_)[compounds],
    MatchStatus = ifelse(compounds == "Unknown EO", "not_found", "resolved"),
    SourceURL = paste0("fixture://cid/", compounds),
    stringsAsFactors = FALSE
  )
  properties = data.frame(
    Query = c("Limonene", "Linalool"),
    CID = c(22311, 6549),
    Title = c("Limonene", "Linalool"),
    MolecularFormula = c("C10H16", "C10H18O"),
    InChIKey = c("XMGQYMWWDOXHJM-UHFFFAOYSA-N",
                 "CDOSHBSSFJOMGT-UHFFFAOYSA-N"),
    CanonicalSMILES = c(NA_character_, "CC(C)=CCCC(C)(C=C)O"),
    IsomericSMILES = c(NA_character_, "CC(C)=CCCC(C)(C=C)O"),
    SMILES = c("CC1=CCC(CC1)C(=C)C", NA_character_),
    ConnectivitySMILES = c("CC1=CCC(CC1)C(=C)C", NA_character_),
    SourceURL = "fixture://properties",
    stringsAsFactors = FALSE
  )
  out = list(
    identity = identity,
    properties = properties,
    synonyms = data.frame(),
    annotations = data.frame(),
    profile = profile
  )
  class(out) = c("uaf_pubchem_profile", class(out))
  out
}

mock_no_smiles_pubchem_profile = function(compounds, ...) {
  list(
    identity = data.frame(Query = compounds,
                          CID = seq_along(compounds),
                          stringsAsFactors = FALSE),
    properties = data.frame(Query = compounds,
                            CID = seq_along(compounds),
                            MolecularFormula = "C",
                            stringsAsFactors = FALSE),
    profile = "ms"
  )
}

write_ainsect_fixture = function(path_id, path_quant) {
  identity = data.frame(
    check.names = FALSE,
    `Melissa ` = c("Linalool", "Limonene", "Linalool"),
    `Clary Sage ` = c("Linalool", "Unknown EO", "")
  )
  quant = data.frame(
    check.names = FALSE,
    Melissa = c(10, 20, 5),
    `Clary Sagae` = c(3, 7, 0)
  )
  utils::write.csv(identity, path_id, row.names = FALSE, na = "")
  utils::write.csv(quant, path_quant, row.names = FALSE, na = "")
}

ainsect_tempdir = function() {
  tmp = tempfile("ainsect-export-")
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  tmp
}

test_that("exportAiNsectMolOlfInputs writes aiNsect CSVs from wide EO tables", {
  tmp = ainsect_tempdir()
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  id_path = file.path(tmp, "ids.csv")
  quant_path = file.path(tmp, "quant.csv")
  out_dir = file.path(tmp, "mololf")
  cache_dir = file.path(tmp, "cache")
  write_ainsect_fixture(id_path, quant_path)

  result = exportAiNsectMolOlfInputs(
    chem_id_csv = id_path,
    chem_quant_csv = quant_path,
    out_dir = out_dir,
    cache_dir = cache_dir,
    profile = "ms",
    throttle = 0,
    uafR_run_id = "test_run",
    profile_fun = mock_ainsect_pubchem_profile
  )

  expect_s3_class(result, "uaf_ainsect_mololf_export")
  expect_true(file.exists(file.path(out_dir, "uafR_compounds.csv")))
  expect_true(file.exists(file.path(out_dir, "treatment_compound_abundance.csv")))
  expect_true(file.exists(file.path(out_dir, "uafR_compounds_unresolved.csv")))
  expect_true(file.exists(file.path(out_dir, "uafR_mololf_export_summary.json")))
  expect_true(file.exists(file.path(out_dir, "uafR_pubchem_identity_audit.csv")))
  expect_true(file.exists(file.path(out_dir, "uafR_pubchem_properties_audit.csv")))

  compounds = utils::read.csv(file.path(out_dir, "uafR_compounds.csv"),
                              stringsAsFactors = FALSE)
  abundance = utils::read.csv(
    file.path(out_dir, "treatment_compound_abundance.csv"),
    stringsAsFactors = FALSE
  )
  unresolved = utils::read.csv(
    file.path(out_dir, "uafR_compounds_unresolved.csv"),
    stringsAsFactors = FALSE
  )

  expect_named(compounds, c("compound_id", "compound_name", "smiles",
                            "inchikey", "pubchem_cid", "molecular_formula",
                            "source", "notes"))
  expect_named(abundance, c("treatment_id", "treatment_name", "compound_id",
                            "relative_abundance", "absolute_abundance",
                            "abundance_units", "uafR_run_id", "source_file"))
  expect_true(all(nzchar(compounds$smiles)))
  expect_true(all(unique(abundance$compound_id) %in% compounds$compound_id))
  expect_equal(sort(unique(abundance$treatment_name)),
               c("Clary Sage", "Melissa"))

  melissa_linalool = abundance[
    abundance$treatment_name == "Melissa" &
      abundance$compound_id == "CDOSHBSSFJOMGT_UHFFFAOYSA_N",
    , drop = FALSE
  ]
  expect_equal(melissa_linalool$absolute_abundance, 15)
  expect_equal(melissa_linalool$relative_abundance, 15 / 35,
               tolerance = 1e-8)

  sums = stats::aggregate(relative_abundance ~ treatment_name,
                          data = abundance, FUN = sum)
  expect_true(all(abs(sums$relative_abundance - 1) < 1e-8))
  expect_equal(unresolved$compound_name, "Unknown EO")
  expect_equal(unresolved$unresolved_reason, "pubchem_query_unresolved")
})

test_that("exportAiNsectMolOlfInputs fails clearly without usable SMILES", {
  tmp = ainsect_tempdir()
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  id_path = file.path(tmp, "ids.csv")
  quant_path = file.path(tmp, "quant.csv")
  write_ainsect_fixture(id_path, quant_path)

  expect_error(
    exportAiNsectMolOlfInputs(
      chem_id_csv = id_path,
      chem_quant_csv = quant_path,
      out_dir = file.path(tmp, "mololf"),
      cache_dir = file.path(tmp, "cache"),
      throttle = 0,
      profile_fun = mock_no_smiles_pubchem_profile
    ),
    "No compounds resolved to usable PubChem SMILES"
  )
})

test_that("exportAiNsectMolOlfInputs rejects unaligned treatments", {
  tmp = ainsect_tempdir()
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  id_path = file.path(tmp, "ids.csv")
  quant_path = file.path(tmp, "quant.csv")
  identity = data.frame(check.names = FALSE, Melissa = "Linalool")
  quant = data.frame(check.names = FALSE, Other = 10)
  utils::write.csv(identity, id_path, row.names = FALSE)
  utils::write.csv(quant, quant_path, row.names = FALSE)

  expect_error(
    exportAiNsectMolOlfInputs(
      chem_id_csv = id_path,
      chem_quant_csv = quant_path,
      out_dir = file.path(tmp, "mololf"),
      cache_dir = file.path(tmp, "cache"),
      throttle = 0,
      profile_fun = mock_ainsect_pubchem_profile
    ),
    "Treatment columns could not be aligned"
  )
})
