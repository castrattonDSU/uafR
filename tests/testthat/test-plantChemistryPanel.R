panel_write_tsv = function(path, header, rows) {
  writeLines(
    c(paste(header, collapse = "\t"),
      vapply(rows, paste, collapse = "\t", FUN.VALUE = character(1))),
    path, useBytes = TRUE
  )
}

panel_resource_fixture = function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  lotus_raw = file.path(root, "lotus_raw.csv")
  lotus_index = file.path(root, "lotus_index.csv")
  utils::write.csv(data.frame(
    traditional_name = "salicin",
    lotus_id = "LTS_PANEL_1",
    allTaxa = "Plantae; Salicaceae; Salix nigra",
    inchikey = "NGFMICBWJRZIBI-UHFFFAOYSA-N",
    smiles = "C1=CC(=C(C=C1CO)O)OC2C(C(C(C(O2)CO)O)O)O",
    cid = "439503",
    stringsAsFactors = FALSE
  ), lotus_raw, row.names = FALSE)
  buildLotusIndex(lotus_raw, out_file = lotus_index, overwrite = TRUE,
                  progress = FALSE)

  raw = file.path(root, "npass_raw")
  dir.create(raw)
  general = file.path(raw, "general.txt")
  structure = file.path(raw, "structure.txt")
  pairs = file.path(raw, "pairs.txt")
  species = file.path(raw, "species.txt")
  panel_write_tsv(
    general,
    c("np_id", "inchikey", "pref_name", "iupac_name", "pubchem_id",
      "molecular_formula"),
    list(c("NPC_PANEL_1", "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "caffeine",
           "caffeine", "2519", "C8H10N4O2"))
  )
  panel_write_tsv(
    structure, c("np_id", "InChI", "InChIKey", "SMILES"),
    list(c("NPC_PANEL_1", "InChI=1S/C8H10N4O2",
           "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
           "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"))
  )
  panel_write_tsv(
    species,
    c("org_id", "org_name", "org_tax_level", "org_tax_id",
      "species_name", "genus_name", "family_name"),
    list(c("ORG_PANEL_1", "Camellia sinensis", "species", "4442",
           "Camellia sinensis", "Camellia", "Theaceae"))
  )
  panel_write_tsv(
    pairs,
    c("src_org_record_id", "src_org_pair_id", "src_org_pair", "org_id",
      "np_id", "new_cp_found", "org_isolation_part",
      "org_collect_location", "org_collect_time", "ref_type", "ref_id",
      "ref_id_type", "ref_url"),
    list(c("REC_PANEL_1", "PAIR_PANEL_1",
           "Camellia sinensis--NPC_PANEL_1", "ORG_PANEL_1", "NPC_PANEL_1",
           "no", "leaf", "China", "2020", "article", "12345678",
           "PMID", "https://pubmed.ncbi.nlm.nih.gov/12345678/"))
  )
  npass_index = file.path(root, "npass_index")
  buildNpassIndex(general, structure, pairs, species, npass_index,
                  progress = FALSE)
  list(lotus = lotus_index, npass = npass_index)
}

panel_release_fixture = function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  version = as.character(utils::packageVersion("uafR"))
  tarball = file.path(root, paste0("uafR_", version, ".tar.gz"))
  writeBin(charToRaw("offline checked source tarball fixture"), tarball)
  manifest = file.path(root, "uafR_release_manifest.json")
  jsonlite::write_json(list(
    workflow = "uafR_release_candidate",
    package = "uafR", package_version = version,
    git_commit = paste(rep("a", 40L), collapse = ""),
    dirty_state_check = "clean",
    source_tarball = basename(tarball),
    source_tarball_path = basename(tarball),
    source_tarball_bytes = as.numeric(file.info(tarball)$size),
    source_tarball_md5 = unname(tools::md5sum(tarball)[[1L]]),
    source_tarball_sha256 = .plant_sha256_file(tarball)
  ), manifest, pretty = TRUE, auto_unbox = TRUE)
  list(manifest = manifest, tarball = tarball)
}

panel_provider_fixtures = function(species) {
  identity_counter = 0L
  identity_keys = c(
    "AAAAAAAAAAAAAA-BBBBBBBBBB-C", "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
    "FFFFFFFFFFFFFF-GGGGGGGGGG-H", "IIIIIIIIIIIIII-JJJJJJJJJJ-K",
    "LLLLLLLLLLLLLL-MMMMMMMMMM-N", "OOOOOOOOOOOOOO-PPPPPPPPPP-Q",
    "RRRRRRRRRRRRRR-SSSSSSSSSS-T", "UUUUUUUUUUUUUU-VVVVVVVVVV-W",
    "XXXXXXXXXXXXXX-YYYYYYYYYY-Z", "BCDEFGHIJKLMNO-CDEFGHIJKL-M",
    "DEFGHIJKLMNOP-EFGHIJKLMN-O", "FGHIJKLMNOPQR-GHIJKLMNOP-Q"
  )
  occurrence = function(provider, compound, prefix, tier =
                         "direct_species_database", confidence = "high",
                         with_identity = TRUE) {
    occurrence_rows = data.frame(
      species = species,
      matched_taxon = species,
      matched_rank = "species",
      compound_name = paste(compound, seq_along(species)),
      source_database = provider,
      source_record_id = paste0(prefix, seq_along(species)),
      evidence_tier = tier,
      confidence = confidence,
      stringsAsFactors = FALSE
    )
    if (!isTRUE(with_identity)) {
      return(list(PlantCompoundOccurrences = occurrence_rows))
    }
    idx = seq.int(identity_counter + 1L,
                  identity_counter + length(species))
    identity_counter <<- max(idx)
    identity_rows = data.frame(
      compound_name = occurrence_rows$compound_name,
      compound_name_clean = .plant_clean_compound(
        occurrence_rows$compound_name
      ),
      source_database = provider,
      source_record_id = occurrence_rows$source_record_id,
      source_compound_id = occurrence_rows$source_record_id,
      source_compound_id_type = paste0(toupper(provider), "_ID"),
      CID = 5000L + idx,
      InChIKey = identity_keys[idx],
      SMILES = paste0(strrep("C", idx), "O"),
      MolecularFormula = paste0("C", idx, "H", 2L * idx + 2L, "O"),
      evidence_url = paste0("https://example.test/", tolower(provider), "/",
                            occurrence_rows$source_record_id),
      evidence_text = "Source-record identity fixture.",
      identity_status = "source_structure_available",
      identity_note = "Offline test fixture; not biological evidence.",
      stringsAsFactors = FALSE
    )
    list(PlantCompoundOccurrences = occurrence_rows,
         SourceCompoundIdentity = identity_rows)
  }
  literature = data.frame(
    species = species,
    source_database = "PubMed",
    source_record_id = paste0("PMID", seq_along(species)),
    pmid = as.character(1000 + seq_along(species)),
    title = paste("Chemical profile fixture for", species),
    evidence_tier = "direct_species_literature",
    confidence = "medium",
    stringsAsFactors = FALSE
  )
  list(
    lotus = occurrence("LOTUS", "lotus compound", "LTS"),
    npass = occurrence("NPASS", "npass compound", "NPS"),
    knapsack = occurrence("KNApSAcK", "knapsack compound", "KNS"),
    pubchem = occurrence("PubChem", "pubchem compound", "PCS"),
    pubmed = list(LiteratureCandidates = literature),
    pubtator = occurrence(
      "PubTator", "pubtator candidate", "PTC",
      tier = "direct_species_pubtator_candidate", confidence = "low",
      with_identity = FALSE
    )
  )
}

panel_enrichment_fixture = function(compounds, detail, query_overrides, ...) {
  cids = seq.int(8001L, length.out = length(compounds))
  list(
    PubChemIdentity = data.frame(
      Query = compounds, CID = cids, MatchStatus = "resolved_fixture",
      stringsAsFactors = FALSE
    ),
    PubChemProperties = data.frame(
      Query = compounds, CID = cids,
      MolecularFormula = paste0("C", seq_along(compounds), "H2O"),
      stringsAsFactors = FALSE
    ),
    PubChemClassifications = data.frame(
      Query = compounds, CID = cids,
      Classification = "terpenoid natural product",
      ClassificationSource = "offline_test_fixture",
      stringsAsFactors = FALSE
    ),
    KEGGMatches = data.frame(
      Query = compounds, KEGG_ID = NA_character_, MatchStatus = "no_records",
      stringsAsFactors = FALSE
    ),
    ValidationSummary = data.frame(
      Status = "pass", ErrorCount = 0L, WarningCount = 0L,
      stringsAsFactors = FALSE
    )
  )
}

panel_server_output_fixture = function(root, species) {
  full = file.path(root, "full")
  dir.create(full, recursive = TRUE, showWarnings = FALSE)
  pairs = utils::combn(sort(species), 2L)
  summary = data.frame(
    group_pair_id = paste(pairs[1L, ], pairs[2L, ], sep = "__"),
    group_a = pairs[1L, ], group_b = pairs[2L, ], group_level = "species",
    compounds_a = 4L, compounds_b = 4L, compound_pair_count = 16L,
    shared_compound_count = 0L, mean_tanimoto = 0.25,
    median_tanimoto = 0.2, p95_tanimoto = 0.7, max_tanimoto = 0.8,
    top_compound_pair_ids = "a | b=0.8",
    top_compound_pairs = "alpha | beta=0.8",
    fingerprint_source = "PubChem_Fingerprint2D",
    compound_pair_count_ge_0_5 = 2L,
    compound_pair_count_ge_0_7 = 1L,
    compound_pair_count_ge_0_85 = 0L,
    compound_pair_count_ge_0_95 = 0L,
    support_status = "computed",
    support_note = "Offline panel fixture.",
    stringsAsFactors = FALSE
  )
  utils::write.csv(summary, file.path(full,
                                      "plant_pair_tanimoto_summary.csv"),
                   row.names = FALSE)
  comparable = function(type, field, value) {
    out = data.frame(
      species_a = pairs[1L, 1L], species_b = pairs[2L, 1L],
      unordered_pair_key = paste(pairs[1L, 1L], pairs[2L, 1L], sep = " || "),
      comparison_type = type, comparison_value = value,
      compound_pair_count = 4L, mean_tanimoto = 0.3,
      median_tanimoto = 0.25, p95_tanimoto = 0.7, max_tanimoto = 0.8,
      compound_pair_count_ge_0_5 = 1L,
      compound_pair_count_ge_0_7 = 1L,
      compound_pair_count_ge_0_85 = 0L,
      compound_pair_count_ge_0_95 = 0L,
      support_tier = "very_low", low_support_caution = "Yes",
      support_note = "Offline panel fixture.",
      comparison_filter = paste0(field, " == '", value, "'"),
      caution = "Source-backed comparable chemistry fixture.",
      stringsAsFactors = FALSE
    )
    out[[field]] = value
    out[, .comparable_pair_summary_cols(type), drop = FALSE]
  }
  utils::write.csv(
    comparable("scope", "comparison_scope", "specialized_metabolites"),
    file.path(full, "comparable_scope_tanimoto_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    comparable("group", "comparison_group", "terpenoids"),
    file.path(full, "comparable_group_tanimoto_summary.csv"),
    row.names = FALSE
  )
  utils::write.csv(
    data.frame(check = "fixture_server_output", status = "pass",
               detail = "Offline test fixture."),
    file.path(full, "server_tanimoto_validation.csv"), row.names = FALSE
  )
  shard = file.path(full, "group_compound_pair_tanimoto.csv.gz")
  con = gzfile(shard, open = "wt")
  utils::write.csv(data.frame(
    group_compound_pair_id = "fixture_pair", group_a = species[[1L]],
    compound_id_a = "a", compound_name_a = "alpha",
    group_b = species[[2L]], compound_id_b = "b", compound_name_b = "beta",
    tanimoto = 0.8, intersection_bits = 8L, union_bits = 10L,
    fingerprint_source = "PubChem_Fingerprint2D"
  ), con, row.names = FALSE)
  close(con)
  utils::write.csv(data.frame(
    mode = "full", run_id = "offline_fixture",
    file = "full/group_compound_pair_tanimoto.csv.gz", row_count = 1L,
    bytes = file.info(shard)$size, md5 = unname(tools::md5sum(shard)),
    stringsAsFactors = FALSE
  ), file.path(root, "server_tanimoto_output_manifest.csv"), row.names = FALSE)
  root
}

test_that("plant chemistry panel completes all-provider offline discovery", {
  root = tempfile("plant_panel_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root)
  species = c("Salix nigra", "Camellia sinensis", "Zea mays")
  input_file = file.path(root, "plants.csv")
  utils::write.csv(data.frame(
    species = species,
    family = c("Salicaceae", "Theaceae", "Poaceae"),
    taxonomy_name_launch_safe = TRUE,
    stringsAsFactors = FALSE
  ), input_file, row.names = FALSE)
  resources = panel_resource_fixture(file.path(root, "resources"))
  release = panel_release_fixture(file.path(root, "release"))
  fixtures = panel_provider_fixtures(species)
  out_dir = file.path(root, "output")
  common = list(
    plant_csv = input_file,
    out_dir = out_dir,
    lotus_index = resources$lotus,
    npass_index = resources$npass,
    sources = c("lotus", "npass", "knapsack", "pubchem", "pubmed",
                "pubtator"),
    taxon_fallback = character(),
    expected_species_count = 3,
    release_manifest = release$manifest,
    source_tarball = release$tarball,
    require_release_artifact = TRUE,
    require_all_providers = TRUE,
    species_chunk_size = 2,
    provider_results = fixtures,
    compound_request_fun = function(url) stop("network forbidden in panel test"),
    progress = FALSE
  )

  preflight = do.call(runPlantChemistryPanel,
                      c(list(mode = "preflight"), common))
  expect_equal(preflight$ExitStatus, 0L)
  indexes = do.call(runPlantChemistryPanel,
                    c(list(mode = "build-indexes"), common))
  expect_equal(indexes$ExitStatus, 0L)
  discovery = do.call(runPlantChemistryPanel,
                      c(list(mode = "discovery"), common))

  expect_equal(discovery$ExitStatus, 0L)
  merged = readRDS(file.path(
    out_dir, "stages", "discovery",
    "plant_phytochemistry_discovery_merged.rds"
  ))
  expect_equal(nrow(merged$PlantQueries), 3L)
  expect_equal(length(unique(merged$PlantQueries$query_id)), 3L)
  expect_setequal(merged$PlantQueries$species, species)
  expect_equal(nrow(merged$ProviderQueryAccounting), 18L)
  expect_equal(length(unique(merged$ProviderQueryAccounting$query_id)), 3L)
  expect_equal(nrow(merged$SourceCompoundIdentity), 12L)
  expect_true(all(merged$ProviderQueryAccounting$query_status %in%
                    c("records", "no_records")))
  expect_true(all(merged$ProviderQueryAccounting$retry_required == "No"))
  gate = utils::read.csv(file.path(
    out_dir, "stages", "discovery", "discovery_quality_gate.csv"
  ), stringsAsFactors = FALSE)
  expect_true(all(gate$status == "pass"))

  writeLines("{\"state\":\"failed\"}", file.path(out_dir, "FAILED.json"))
  writeLines("{\"state\":\"paused\"}",
             file.path(out_dir, "PAUSED_SERVICE_BUSY.json"))
  reused = do.call(runPlantChemistryPanel,
                   c(list(mode = "discovery"), common))
  expect_equal(reused$StageResults$discovery$status, "reused")
  expect_false(file.exists(file.path(out_dir, "FAILED.json")))
  expect_false(file.exists(file.path(out_dir, "PAUSED_SERVICE_BUSY.json")))

  pilot = do.call(runPlantChemistryPanel, c(
    list(mode = "pilot", enrichment_fun = panel_enrichment_fixture), common
  ))
  expect_equal(pilot$ExitStatus, 0L)
  pilot_result = readRDS(file.path(
    out_dir, "stages", "pilot", "pilot_result.rds"
  ))
  expect_gt(nrow(pilot_result$CompoundResolution), 0L)
  expect_true(all(pilot_result$CompoundResolution$resolution_source ==
                    "not_attempted"))
  expect_false(any(pilot_result$CompoundResolution$resolved %in% TRUE))
  pilot_selection = utils::read.csv(file.path(
    out_dir, "stages", "pilot", "pilot_enrichment_selection.csv"
  ), stringsAsFactors = FALSE)
  expect_gt(nrow(pilot_selection), 0L)
  expect_lte(nrow(pilot_selection), 25L)
  expect_true(all(pilot_selection$SourceIdentityAvailable == "Yes"))
  expect_true(all(pilot_selection$QueryType %in%
                    c("source_cid", "source_inchikey")))
  pilot_gate = utils::read.csv(file.path(
    out_dir, "stages", "pilot", "pilot_quality_gate.csv"
  ), stringsAsFactors = FALSE)
  expect_true(all(pilot_gate$status == "pass"))
  expect_true(all(c("pilot_discovery_record_limit",
                    "pilot_enrichment_source_identity") %in%
                  pilot_gate$check))
  identity = do.call(runPlantChemistryPanel,
                     c(list(mode = "identity"), common))
  expect_equal(identity$ExitStatus, 0L)
  research = do.call(runPlantChemistryPanel, c(
    list(mode = "research-enrichment",
         enrichment_fun = panel_enrichment_fixture), common
  ))
  expect_equal(research$ExitStatus, 0L)
  full = do.call(runPlantChemistryPanel, c(
    list(mode = "full-enrichment",
         enrichment_fun = panel_enrichment_fixture), common
  ))
  expect_equal(full$ExitStatus, 0L)
  handoff = do.call(runPlantChemistryPanel,
                    c(list(mode = "tanimoto-handoff"), common))
  expect_equal(handoff$ExitStatus, 0L)
  expect_true(file.exists(file.path(
    out_dir, "tanimoto_handoff", "plant_species_universe.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir, "tanimoto_handoff", "run_plant_tanimoto_server.R"
  )))
  expect_true(file.exists(file.path(
    out_dir, "tanimoto_handoff", basename(release$manifest)
  )))
  expect_true(file.exists(file.path(
    out_dir, "tanimoto_handoff", basename(release$tarball)
  )))
  command_text = readLines(file.path(
    out_dir, "tanimoto_handoff", "RUN_TANIMOTO_SERVER.txt"
  ), warn = FALSE)
  expect_true(any(grepl("--species-universe", command_text, fixed = TRUE)))
  expect_true(any(grepl("--release-manifest", command_text, fixed = TRUE)))
  expect_true(any(grepl("R CMD INSTALL", command_text, fixed = TRUE)))
  source_manifest = utils::read.csv(file.path(
    out_dir, "tanimoto_handoff", "source_release_manifest.csv"
  ), stringsAsFactors = FALSE)
  expect_true(all(!grepl("^/", source_manifest$source_path)))
  expect_true(all(c("release_manifest", "source_tarball") %in%
                    source_manifest$role))

  server_root = panel_server_output_fixture(
    file.path(out_dir, "tanimoto_server_output"), species
  )
  finalized = do.call(runPlantChemistryPanel, c(
    list(mode = "finalize", server_results_dir = server_root), common
  ))
  expect_equal(finalized$ExitStatus, 0L)
  reconciliation = utils::read.csv(file.path(
    out_dir, "plant_chemistry_analysis_bundle",
    "99b_FinalReconciliation.csv"
  ), stringsAsFactors = FALSE)
  expect_true(all(reconciliation$status == "pass"))
  expect_true(file.exists(file.path(
    out_dir, "plant_chemistry_analysis_bundle",
    "22_ComparableScopeTanimotoSummary.csv"
  )))
})

test_that("panel pilot limits discovery and selects source identities", {
  expect_equal(.plant_panel_pilot_record_limit(Inf, 25L), 25L)
  expect_equal(.plant_panel_pilot_record_limit(10L, 25L), 10L)
  expect_equal(.plant_panel_pilot_record_limit(Inf, 40L), 40L)

  source = data.frame(
    compound_name = c("Beta", "Alpha", "Gamma", "Alpha"),
    CID = c(22L, NA, 33L, NA),
    InChIKey = c(
      NA, "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
      "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
      "AAAAAAAAAAAAAA-BBBBBBBBBB-C"
    ),
    SMILES = c("CC", "CO", "CCC", "CO"),
    stringsAsFactors = FALSE
  )
  selected = .plant_panel_pilot_enrichment_input(source, 2L)
  expect_equal(nrow(selected), 2L)
  expect_equal(selected$QueryType, c("source_cid", "source_cid"))
  expect_equal(selected$compound_name, c("Beta", "Gamma"))
  expect_true(all(selected$SourceIdentityAvailable == "Yes"))
  expect_equal(length(unique(selected$PubChemQuery)), nrow(selected))
})

test_that("research enrichment selection is bounded and species aware", {
  occurrences = .plant_normalize_occurrences(data.frame(
    species = c("Plant one", "Plant two", "Plant three", "Plant one"),
    matched_taxon = c("Plant one", "Plant two", "Plant three", "Plant"),
    matched_rank = c("species", "species", "species", "genus"),
    compound_name = c("caffeine", "salicin", "limonene", "quercetin"),
    source_database = c("NPASS", "LOTUS", "LOTUS", "LOTUS"),
    source_record_id = c("N1", "L1", "L2", "L3"),
    evidence_tier = c(
      "direct_species_database", "direct_species_database",
      "direct_species_database", "genus_database_fallback"
    ),
    confidence = c("high", "high", "high", "medium"),
    stringsAsFactors = FALSE
  ))
  resolution = data.frame(
    compound_name = c("caffeine", "salicin", "limonene", "quercetin"),
    compound_name_clean = c("caffeine", "salicin", "limonene", "quercetin"),
    query_count = 1L,
    resolved = TRUE,
    CID = c("2519", "439503", "22311", "5280343"),
    InChIKey = NA_character_,
    SMILES = c("CN", "CO", "CC", "C1=CC=CC=C1"),
    MolecularFormula = c("C8H10N4O2", "C13H18O7", "C10H16", "C15H10O7"),
    resolution_source = "source_record",
    notes = NA_character_,
    stringsAsFactors = FALSE
  )
  result = list(
    CompoundResolution = resolution,
    PlantCompoundOccurrences = occurrences
  )

  first = .plant_panel_enrichment_selection(result, limit = 3L)
  second = .plant_panel_enrichment_selection(result, limit = 3L)

  expect_equal(first$included$compound_name_clean,
               second$included$compound_name_clean)
  expect_equal(nrow(first$included), 3L)
  expect_setequal(first$included$compound_name_clean,
                  c("caffeine", "salicin", "limonene"))
  expect_true(all(first$included$research_priority_reason ==
                    "species_coverage"))
  expect_equal(first$Summary$eligible_species_count, 3L)
  expect_equal(first$Summary$selected_species_count, 3L)
  expect_equal(first$Summary$selected_species_coverage_fraction, 1)
  expect_equal(first$Summary$research_deferred_eligible_rows, 1L)
  expect_equal(
    first$excluded$enrichment_exclusion_reason[
      first$excluded$compound_name_clean == "quercetin"
    ],
    "research_enrichment_limit"
  )

  expanded = .plant_panel_enrichment_selection(result, limit = 4L)
  expect_equal(expanded$included$research_priority_reason[[4L]],
               "global_evidence_priority")
  expect_equal(expanded$included$compound_name_clean[[4L]], "quercetin")
})

test_that("pilot resume clears derived exports but preserves checkpoints", {
  root = tempfile("pilot_resume_exports_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  checkpoint_dir = file.path(root, "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE)
  writeLines("derived", file.path(root, "all_occurrences.csv"))
  writeLines("checkpoint", file.path(checkpoint_dir, "chunk.rds"))

  removed = .plant_panel_clear_derived_exports(root)

  expect_true(file.path(root, "all_occurrences.csv") %in% removed)
  expect_false(file.exists(file.path(root, "all_occurrences.csv")))
  expect_true(file.exists(file.path(checkpoint_dir, "chunk.rds")))
})

test_that("fallback selection includes only provider-specific direct no-hits", {
  result = list(ProviderQueryAccounting = data.frame(
    provider = c("lotus", "lotus", "npass"),
    query_status = c("records", "no_records", "no_records"),
    species = c("Plant one", "Plant two", "Plant three"),
    query_plant = c("Plant one", "Plant two", "Plant three"),
    stringsAsFactors = FALSE
  ))

  expect_equal(.plant_panel_provider_no_hit_species(result, "lotus"),
               "Plant two")
  expect_equal(.plant_panel_provider_no_hit_species(result, "npass"),
               "Plant three")
  expect_length(.plant_panel_provider_no_hit_species(NULL, "lotus"), 0L)
})

test_that("panel markers reject changed supporting inputs", {
  root = tempfile("plant_panel_signature_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root)
  plant_file = file.path(root, "plants.csv")
  support_file = file.path(root, "review.csv")
  exclusion_file = file.path(root, "excluded.csv")
  utils::write.csv(data.frame(species = c("Plant one", "Plant two")),
                   plant_file, row.names = FALSE)
  utils::write.csv(data.frame(note = "first"), support_file,
                   row.names = FALSE)
  utils::write.csv(data.frame(species = "Plant excluded"), exclusion_file,
                   row.names = FALSE)
  first = runPlantChemistryPanel(
    "preflight", plant_file, file.path(root, "output"),
    exclusion_files = exclusion_file,
    supporting_files = support_file, sources = "pubmed",
    require_all_providers = FALSE, progress = FALSE
  )
  expect_equal(first$ExitStatus, 0L)
  manifest = utils::read.csv(file.path(
    root, "output", "inputs", "input_manifest.csv"
  ), stringsAsFactors = FALSE)
  expect_equal(manifest$role[match(basename(plant_file),
                                   manifest$snapshot_file)], "plant_input")
  expect_equal(manifest$role[match(basename(exclusion_file),
                                   manifest$snapshot_file)],
               "exclusion_ledger")
  expect_equal(manifest$role[match(basename(support_file),
                                   manifest$snapshot_file)],
               "supporting_ledger")
  utils::write.csv(data.frame(note = "changed"), support_file,
                   row.names = FALSE)
  expect_error(
    runPlantChemistryPanel(
      "preflight", plant_file, file.path(root, "output"),
      exclusion_files = exclusion_file,
      supporting_files = support_file, sources = "pubmed",
      require_all_providers = FALSE, progress = FALSE
    ),
    "Input snapshot differs"
  )
})

test_that("production release artifacts are commit and checksum bound", {
  root = tempfile("panel_release_contract_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  release = panel_release_fixture(root)
  config = list(
    release_manifest = release$manifest,
    source_tarball = release$tarball,
    require_release_artifact = TRUE,
    package_version = as.character(utils::packageVersion("uafR"))
  )
  valid = .plant_panel_validate_release_artifacts(config)
  expect_true(all(valid$Checks$status == "pass"))
  expect_match(valid$GitCommit, "^[a-f0-9]{40}$")

  writeBin(charToRaw("tampered source tarball"), release$tarball)
  invalid = .plant_panel_validate_release_artifacts(config)
  expect_equal(invalid$Checks$status[
    invalid$Checks$check == "release_manifest_tarball_md5"
  ], "fail")
  expect_equal(invalid$Checks$status[
    invalid$Checks$check == "release_manifest_tarball_sha256"
  ], "fail")

  missing = .plant_panel_validate_release_artifacts(list(
    release_manifest = NA_character_, source_tarball = NA_character_,
    require_release_artifact = TRUE,
    package_version = as.character(utils::packageVersion("uafR"))
  ))
  expect_equal(missing$Checks$status, "fail")
})

test_that("NPASS production download URLs retain source-role names", {
  urls = uafR:::.plant_panel_npass_urls()

  expect_named(
    urls,
    c("general_info", "structures", "species_pairs", "species_info")
  )
  expect_match(urls[["general_info"]],
               "NPASS3[.]0_naturalproducts_generalinfo[.]txt$")
  expect_match(urls[["structures"]],
               "NPASS3[.]0_naturalproducts_structure[.]txt$")
  expect_match(urls[["species_pairs"]],
               "NPASS3[.]0_naturalproducts_species_pair[.]txt$")
  expect_match(urls[["species_info"]], "NPASS3[.]0_species_info[.]txt$")
})

test_that("production resource downloads resume persistent partial files", {
  root = tempfile("panel_download_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root)
  destination = file.path(root, "resource.txt")
  offsets = numeric()
  attempts = 0L
  transfer = function(url, partial, offset, request_timeout, expected) {
    offsets <<- c(offsets, offset)
    attempts <<- attempts + 1L
    connection = file(partial, open = if (offset > 0) "ab" else "wb")
    on.exit(close(connection), add = TRUE)
    if (attempts == 1L) {
      writeBin(charToRaw("1234"), connection)
      stop("simulated interrupted transfer")
    }
    writeBin(charToRaw("567890"), connection)
  }

  uafR:::.plant_panel_download_resource(
    "https://example.org/resource.txt", destination,
    request_timeout = 60, max_attempts = 2, retry_wait = 0,
    progress = FALSE,
    metadata_fun = function(url, timeout) list(expected_bytes = 10),
    transfer_fun = transfer, sleep_fun = function(seconds) NULL
  )

  expect_equal(offsets, c(0, 4))
  expect_equal(readChar(destination, nchars = 10, useBytes = TRUE),
               "1234567890")
  expect_false(file.exists(paste0(destination, ".partial")))

  uafR:::.plant_panel_download_resource(
    "https://example.org/resource.txt", destination,
    request_timeout = 60, max_attempts = 1, retry_wait = 0,
    progress = FALSE,
    metadata_fun = function(url, timeout) list(expected_bytes = 10),
    transfer_fun = function(...) stop("complete files must be reused"),
    sleep_fun = function(seconds) NULL
  )
})

test_that("panel distinguishes completed, service-busy, and failed discovery", {
  make_result = function(run_status, discovery_complete, chunk_status,
                         error_message = NA_character_) {
    list(
      BatchRunManifest = data.frame(
        run_status = run_status,
        discovery_complete = discovery_complete,
        pause_reason = error_message,
        stringsAsFactors = FALSE
      ),
      BatchChunkManifest = data.frame(
        status = chunk_status,
        error_message = error_message,
        stringsAsFactors = FALSE
      )
    )
  }

  completed = .plant_panel_batch_completion_state(
    make_result("completed", "Yes", "completed")
  )
  paused = .plant_panel_batch_completion_state(
    make_result("incomplete", "No", "rate_limited",
                "HTTP 500 temporarily unavailable")
  )
  failed = .plant_panel_batch_completion_state(
    make_result("incomplete", "No", "failed", "malformed response")
  )

  expect_equal(completed$state, "completed")
  expect_equal(paused$state, "paused_service_busy")
  expect_equal(failed$state, "incomplete")
})
