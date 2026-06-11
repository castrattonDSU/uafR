plant_enrichment_fixture = function(compounds, ...) {
  data.frame_for = function(name, cid, key, smiles, formula) {
    data.frame(Query = name,
               CID = cid,
               MolecularFormula = formula,
               InChIKey = key,
               CanonicalSMILES = smiles,
               IsomericSMILES = smiles,
               stringsAsFactors = FALSE)
  }
  props = do.call(rbind, lapply(seq_along(compounds), function(i) {
    compound = compounds[[i]]
    data.frame_for(compound,
                   1000 + i,
                   paste0("TESTINCHIKEY", i),
                   paste0("C", i, "CO"),
                   paste0("C", i + 5, "H", i + 10, "O"))
  }))
  traits = data.frame(
    Query = compounds,
    TraitType = "chemistry",
    TraitGroup = c("phenolic", rep("volatile_proxy", length(compounds) - 1)),
    TraitValue = c("phenolic glycoside", rep("volatile candidate",
                                             length(compounds) - 1)),
    Confidence = "high",
    SourceDatabase = "fixture",
    stringsAsFactors = FALSE
  )
  evidence = data.frame(
    Query = compounds,
    CID = props$CID,
    EvidenceType = "trait",
    AnalysisKey = paste0("trait__", seq_along(compounds)),
    SourceDatabase = "fixture",
    EvidenceText = paste("Fixture evidence for", compounds),
    EvidenceURL = "https://example.test/evidence",
    Confidence = "high",
    stringsAsFactors = FALSE
  )
  classes = data.frame(
    Query = compounds,
    ClassName = c("Phenolic glycosides", rep("Volatile organic compounds",
                                             length(compounds) - 1)),
    ClassGroup = c("phenolics", rep("volatiles", length(compounds) - 1)),
    Superclass = "Organic compounds",
    Class = c("Phenolics", rep("Volatiles", length(compounds) - 1)),
    Subclass = c("Glycosides", rep("Terpenoid-like", length(compounds) - 1)),
    stringsAsFactors = FALSE
  )
  list(PubChemProperties = props,
       ChemicalTraits = traits,
       ChemicalTraitEvidence = evidence,
       ChemicalClasses = classes,
       ValidationSummary = data.frame(Status = "pass", stringsAsFactors = FALSE))
}

plant_pubchem_profile_fixture = function(compounds, profile = "safety", ...) {
  out = .categorate_empty_pubchem_profile(compounds, profile)
  cids = seq_along(compounds) + 9000L
  out$identity = data.frame(
    Query = compounds,
    CID = cids,
    MatchStatus = "resolved",
    SourceURL = paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cids),
    stringsAsFactors = FALSE
  )
  out$properties = data.frame(
    Query = compounds,
    CID = cids,
    Title = compounds,
    MolecularFormula = c("C13H18O7", "C8H10N4O2")[seq_along(compounds)],
    MolecularWeight = c(286.28, 194.19)[seq_along(compounds)],
    IUPACName = compounds,
    InChI = paste0("InChI=1S/fixture", seq_along(compounds)),
    InChIKey = paste0("PUBCHEMFIXTURE", seq_along(compounds)),
    CanonicalSMILES = c("C1=CC=C(C=C1)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")[seq_along(compounds)],
    IsomericSMILES = c("C1=CC=C(C=C1)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")[seq_along(compounds)],
    SMILES = c("C1=CC=C(C=C1)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")[seq_along(compounds)],
    ConnectivitySMILES = c("C1=CC=C(C=C1)O", "Cn1cnc2c1c(=O)n(C)c(=O)n2C")[seq_along(compounds)],
    ExactMass = c(286.105, 194.080)[seq_along(compounds)],
    MonoisotopicMass = c(286.105, 194.080)[seq_along(compounds)],
    XLogP = c(0.4, -0.1)[seq_along(compounds)],
    TPSA = c(120, 61)[seq_along(compounds)],
    Complexity = c(300, 293)[seq_along(compounds)],
    Charge = 0,
    HBondDonorCount = c(4, 0)[seq_along(compounds)],
    HBondAcceptorCount = c(7, 6)[seq_along(compounds)],
    RotatableBondCount = c(3, 0)[seq_along(compounds)],
    HeavyAtomCount = c(20, 14)[seq_along(compounds)],
    IsotopeAtomCount = 0,
    AtomStereoCount = 0,
    DefinedAtomStereoCount = 0,
    UndefinedAtomStereoCount = 0,
    BondStereoCount = 0,
    DefinedBondStereoCount = 0,
    UndefinedBondStereoCount = 0,
    CovalentUnitCount = 1,
    SourceURL = paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", cids),
    stringsAsFactors = FALSE
  )
  out$provenance = data.frame(
    Table = "properties",
    SourceURL = "https://pubchem.ncbi.nlm.nih.gov/",
    stringsAsFactors = FALSE
  )
  out
}

plant_comparability_enrichment_fixture = function(compounds, ...) {
  props = data.frame(
    Query = compounds,
    CID = seq_along(compounds) + 7000L,
    MolecularFormula = paste0("C", seq_along(compounds) + 5, "H12O2"),
    InChIKey = paste0("COMPARABILITYKEY", seq_along(compounds)),
    CanonicalSMILES = paste0("C", seq_along(compounds), "O"),
    IsomericSMILES = paste0("C", seq_along(compounds), "O"),
    stringsAsFactors = FALSE
  )
  trait_map = data.frame(
    Query = c("salicin", "caffeine", "limonene", "glucose",
              "linoleic acid", "alliin", "abscisic acid"),
    TraitType = "chemistry",
    TraitGroup = c("phenolic", "alkaloid", "terpenoid", "carbohydrate",
                   "fatty acid", "organosulfur", "plant hormone"),
    TraitValue = c("phenolic glycoside", "alkaloid nitrogenous compound",
                   "volatile monoterpene", "primary carbohydrate",
                   "fatty acid lipid", "organosulfur sulfoxide",
                   "plant hormone signaling compound"),
    Confidence = "high",
    SourceDatabase = "fixture",
    stringsAsFactors = FALSE
  )
  trait_map = trait_map[trait_map$Query %in% compounds, , drop = FALSE]
  class_map = data.frame(
    Query = trait_map$Query,
    ClassName = c("Phenolic glycosides", "Alkaloids", "Monoterpenes",
                  "Carbohydrates", "Fatty acids", "Organosulfur compounds",
                  "Plant hormones")[match(
                    trait_map$Query,
                    c("salicin", "caffeine", "limonene", "glucose",
                      "linoleic acid", "alliin", "abscisic acid")
                  )],
    ClassGroup = trait_map$TraitGroup,
    Superclass = "Organic compounds",
    Class = trait_map$TraitGroup,
    Subclass = trait_map$TraitValue,
    stringsAsFactors = FALSE
  )
  evidence = data.frame(
    Query = trait_map$Query,
    CID = seq_len(nrow(trait_map)) + 7000L,
    EvidenceType = "trait",
    AnalysisKey = paste0("comparison_trait__", seq_len(nrow(trait_map))),
    SourceDatabase = "fixture",
    EvidenceText = paste("Fixture evidence for", trait_map$Query),
    EvidenceURL = "https://example.test/comparability",
    Confidence = "high",
    stringsAsFactors = FALSE
  )
  list(PubChemProperties = props,
       ChemicalTraits = trait_map,
       ChemicalTraitEvidence = evidence,
       ChemicalClasses = class_map,
       ValidationSummary = data.frame(Status = "pass", stringsAsFactors = FALSE))
}

plant_rich_class_enrichment_fixture = function(compounds, ...) {
  props = data.frame(
    Query = compounds,
    CID = seq_along(compounds) + 8000L,
    MolecularFormula = paste0("C", seq_along(compounds) + 10, "H16O3"),
    InChIKey = paste0("RICHCLASSKEY", seq_along(compounds)),
    CanonicalSMILES = paste0("CC", seq_along(compounds), "O"),
    IsomericSMILES = paste0("CC", seq_along(compounds), "O"),
    stringsAsFactors = FALSE
  )
  classes = data.frame(
    Query = c("plain phenolic", "glucose"),
    CID = c(8001L, 8006L),
    ClassSystem = "LOTUS",
    ClassType = "natural_product_class",
    ClassID = c("LTS_CLASS_001", "LTS_CLASS_002"),
    ClassName = c("Shikimates and Phenylpropanoids", "Alkaloids"),
    ClassGroup = c("shikimates phenylpropanoids", "alkaloids"),
    ClassPath = c("Organic compounds > Shikimates and Phenylpropanoids",
                  "Organic compounds > Alkaloids"),
    SourceTable = "ChemicalClasses",
    EvidenceText = c("LOTUS class tree reports Shikimates and Phenylpropanoids",
                     "Contrived source-backed class row for priority test"),
    EvidenceURL = "https://example.test/classes",
    ExtractionRule = "fixture_source_class",
    Confidence = "high",
    stringsAsFactors = FALSE
  )
  lotus = data.frame(
    Query = "plain alkaloid",
    CID = 8003L,
    LOTUS_ID = "LTS000003",
    NaturalProductClass = "Alkaloids",
    RawValue = "Alkaloids reported by LOTUS profile",
    SourceURL = "https://example.test/lotus",
    stringsAsFactors = FALSE
  )
  pubchem_classes = data.frame(
    Query = "plain terpene",
    CID = 8002L,
    Source = "LOTUS",
    TreeID = "123",
    TreeName = "LOTUS natural products",
    TreeType = "chemical",
    RootHNID = NA_character_,
    HNID = "456",
    NodeID = "456",
    ParentNodeID = "455",
    ClassName = "Monoterpenoids",
    ParentClass = "Prenol lipids",
    ClassPath = "Organic compounds > Prenol lipids > Monoterpenoids",
    ClassDepth = 3,
    CompoundCount = NA_integer_,
    TaxonomyCount = NA_integer_,
    DOICount = NA_integer_,
    PubMedCount = NA_integer_,
    SourceURL = "https://example.test/pubchem/source",
    ClassificationURL = "https://example.test/pubchem/class",
    PubChemURL = "https://pubchem.ncbi.nlm.nih.gov/compound/8002",
    stringsAsFactors = FALSE
  )
  ontology = data.frame(
    Query = "plain polyketide",
    CID = 8005L,
    OntologyDomain = "ecology",
    OntologyGroup = "natural_product_class",
    OntologyTerm = "polyketides",
    OntologyLabel = "Polyketides",
    OntologyKey = "ecology__natural_product_class__polyketides",
    OntologyID = NA_character_,
    OntologyIDSource = NA_character_,
    OntologyIDURL = NA_character_,
    SourceTraitType = "chemical_class",
    SourceTraitGroup = "natural_product_class",
    SourceTraitValue = "Polyketides",
    SourceTraitKey = "chemical_class__natural_product_class__polyketides",
    SourceDatabase = "fixture",
    Confidence = "high",
    ConfidenceScore = 1,
    stringsAsFactors = FALSE
  )
  kegg = data.frame(
    Query = "plain lipid",
    KEGG_ID = "C00001",
    Database = "compound",
    ClassificationType = "pathway_group",
    Classification = "Lipid metabolism; Fatty acid biosynthesis",
    Evidence = "KEGG pathway group fixture",
    EvidenceURL = "https://www.kegg.jp/pathway/map01212",
    RetrievedAt = "2026-01-01T00:00:00+0000",
    stringsAsFactors = FALSE
  )
  derived = data.frame(
    Query = "plain volatile",
    CID = 8007L,
    natural_product_superclasses = NA_character_,
    natural_product_classes = "Terpenoids",
    natural_product_subclasses = "Monoterpenoids",
    metabolic_context = NA_character_,
    ecological_context = "Terpenoids",
    analytical_context = "high volatility proxy",
    volatility_proxy = "high",
    lipophilicity_bin = "balanced",
    polarity_bin = "low",
    stringsAsFactors = FALSE
  )
  list(
    PubChemProperties = props,
    ChemicalClasses = classes,
    LOTUSProfile = lotus,
    PubChemClassifications = pubchem_classes,
    ChemicalTraitOntology = ontology,
    KEGGClassifications = kegg,
    DerivedGroups = derived,
    ChemicalTraits = .uaf_empty_table(c("Query", "CID", "TraitType",
                                        "TraitGroup", "TraitValue",
                                        "Confidence", "SourceDatabase")),
    ChemicalTraitEvidence = .uaf_empty_table(c("Query", "CID",
                                               "EvidenceType",
                                               "AnalysisKey",
                                               "SourceDatabase",
                                               "EvidenceText",
                                               "EvidenceURL",
                                               "Confidence")),
    ValidationSummary = data.frame(Status = "pass", stringsAsFactors = FALSE)
  )
}

pubmed_fixture_request = function(url) {
  if (grepl("esearch", url)) {
    return(list(esearchresult = list(idlist = list("111", "222"))))
  }
  if (grepl("esummary", url)) {
    return(list(result = list(
      `111` = list(title = "Phytochemical profile of Salix nigra bark",
                   articleids = list(list(idtype = "doi",
                                          value = "10.1000/salix"))),
      `222` = list(title = "Metabolites reported from Salix nigra roots",
                   articleids = list())
    )))
  }
  list()
}

pubtator_fixture_request = function(url) {
  list(documents = list(list(
    pmid = "333",
    title = "Salix nigra phytochemical annotations",
    abstract = "Salicin was discussed with Salix nigra in a chemical ecology context.",
    annotations = list(
      list(text = "Salix nigra", infons = list(type = "Species")),
      list(text = "salicin", infons = list(type = "Chemical"))
    )
  )))
}

test_that("plantPhytochemistrySchema returns core table contracts", {
  schema = plantPhytochemistrySchema()
  expect_true(all(c("PlantQueries", "PlantCompoundOccurrences",
                    "PlantContextEvidence",
                    "SpeciesChemistrySummary", "ChemistryComparability",
                    "ComparableChemistryMatrix", "Provenance") %in%
                    schema$Table))
  occurrence = plantPhytochemistrySchema("PlantCompoundOccurrences")
  expect_true(all(c("species", "compound_name", "evidence_tier",
                    "confidence", "occurrence_status", "analysis_ready",
                    "plant_part_group", "evidence_quality_score") %in%
                    occurrence$Column))
  context = plantPhytochemistrySchema("PlantContextEvidence")
  expect_true(all(c("context_type", "raw_context_text",
                    "normalized_context", "source_field",
                    "extraction_rule", "context_confidence") %in%
                    context$Column))
  comparability = plantPhytochemistrySchema("ChemistryComparability")
  expect_true(all(c("metabolism_domain", "biosynthetic_family",
                    "chemical_behavior", "comparison_scope",
                    "comparability_confidence",
                    "comparable_for_matrix") %in%
                    comparability$Column))
})

test_that("plant name parsing canonicalizes query case and flags ambiguous taxa", {
  phyto = resolvePlantPhytochemistry(
    plants = c("Acacia Erioloba", "Artocarpus Sp.", "Acer",
               "Vitis Berlandieri X Riparia"),
    sources = character(),
    enrich_compounds = FALSE,
    detail = "none"
  )

  expect_equal(phyto$PlantQueries$query_plant[[1]], "Acacia Erioloba")
  expect_equal(phyto$PlantQueries$species[[1]], "Acacia erioloba")
  expect_equal(phyto$PlantQueries$species[[4]], "Vitis berlandieri x riparia")
  expect_equal(
    phyto$PlantNameResolution$query_status,
    c("parsed_species", "name_needs_review", "name_needs_review",
      "parsed_species")
  )
  expect_match(phyto$PlantNameResolution$resolution_note[[1]],
               "Species-like binomial")
})

test_that("LOTUS prefix reader returns bounded valid JSON", {
  txt = paste0(
    '{"originalQuery":"Achillea millefolium",',
    '"determinedInputType":"name","naturalProducts":[',
    '{"lotus_id":"LTS1","traditional_name":"alpha",',
    '"taxonomyReferenceObjects":{"r1":{"NCBI":[{"species":"Achillea millefolium"}]}}},',
    '{"lotus_id":"LTS2","traditional_name":"beta with } in text",',
    '"fragments":{"[C]({nested})":1},',
    '"taxonomyReferenceObjects":{"r2":{"NCBI":[{"species":"Achillea millefolium"}]}}},',
    '{"lotus_id":"LTS3","traditional_name":"gamma",',
    '"taxonomyReferenceObjects":{"r3":{"NCBI":[{"species":"Achillea millefolium"}]}}}',
    ']}'
  )
  prefix = .plant_lotus_json_prefix(txt, max_records = 2)
  parsed = jsonlite::fromJSON(prefix, simplifyVector = FALSE)
  expect_equal(length(parsed$naturalProducts), 2)
  expect_equal(parsed$naturalProducts[[1]]$lotus_id, "LTS1")
  expect_equal(parsed$naturalProducts[[2]]$lotus_id, "LTS2")
})

test_that("LOTUS live record cap is conservative and configurable", {
  old = Sys.getenv("UAFR_LOTUS_MAX_RECORDS_PER_SPECIES", unset = NA)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("UAFR_LOTUS_MAX_RECORDS_PER_SPECIES")
    } else {
      Sys.setenv(UAFR_LOTUS_MAX_RECORDS_PER_SPECIES = old)
    }
  }, add = TRUE)

  Sys.unsetenv("UAFR_LOTUS_MAX_RECORDS_PER_SPECIES")
  expect_equal(.plant_lotus_live_record_cap(100), 2)
  expect_equal(.plant_lotus_live_record_cap(1), 1)

  Sys.setenv(UAFR_LOTUS_MAX_RECORDS_PER_SPECIES = "4")
  expect_equal(.plant_lotus_live_record_cap(100), 4)
  expect_match(.plant_lotus_live_cap_note(100, 4, TRUE),
               "capped at 4 record")
})

test_that("local LOTUS indexes expand taxon fields and label fallback evidence", {
  lotus = data.frame(
    traditional_name = c("Salicin", "Salixalbin", "Achilleol A"),
    lotus_id = c("LTS000010", "LTS000012", "LTS000011"),
    allTaxa = c("Plantae; Salicaceae; Salix nigra | Salix alba",
                "Plantae; Salicaceae; Salix alba",
                "Plantae; Asteraceae; Achillea millefolium"),
    chemicalTaxonomyNPclassifierPathway = c("Shikimates and Phenylpropanoids",
                                            "Terpenoids",
                                            "Terpenoids"),
    chemicalTaxonomyNPclassifierSuperclass = c("Phenolic glycosides",
                                               "Diterpenoids",
                                               "Sesquiterpenoids"),
    doi = c("10.1000/salix", "10.1000/alba", "10.1000/achillea"),
    stringsAsFactors = FALSE
  )

  index = standardizeLotusIndex(lotus)
  tmp_index = tempfile(fileext = ".csv")
  utils::write.csv(lotus, tmp_index, row.names = FALSE)
  expect_equal(nrow(standardizeLotusIndex(tmp_index)), nrow(index))
  expect_true(all(c("Salix nigra", "Salix alba",
                    "Achillea millefolium") %in% index$species))
  expect_equal(index$source_database[[1]], "LOTUS")
  expect_true(any(index$chemical_class_pathway ==
                    "Shikimates and Phenylpropanoids"))

  occ = queryLotusIndex(
    data.frame(species = "Salix nigra", family = "Salicaceae",
               stringsAsFactors = FALSE),
    index,
    taxon_fallback = c("species", "genus")
  )
  expect_true(any(occ$compound_name == "Salicin" &
                    occ$matched_rank == "species"))
  expect_true(any(occ$compound_name == "Salixalbin" &
                    occ$matched_rank == "genus"))
  expect_true(any(occ$analysis_ready[occ$matched_rank == "species"] == "Yes"))
  expect_true(all(occ$analysis_ready[occ$matched_rank == "genus"] == "No"))
})

test_that("resolver uses a local LOTUS index without live requests", {
  lotus = data.frame(
    traditional_name = "Achilleol A",
    lotus_id = "LTS000011",
    allTaxa = "Plantae; Asteraceae; Achillea millefolium",
    stringsAsFactors = FALSE
  )
  called = FALSE
  phyto = resolvePlantPhytochemistry(
    "Achillea millefolium",
    sources = "lotus",
    lotus_index = lotus,
    enrich_compounds = FALSE,
    detail = "none",
    request_fun = function(url) {
      called <<- TRUE
      stop("live LOTUS request should not be used")
    }
  )

  expect_false(called)
  expect_equal(nrow(phyto$PlantCompoundOccurrences), 1)
  expect_equal(phyto$PlantCompoundOccurrences$source_database, "LOTUS")
  expect_equal(phyto$ProviderDiagnostics$request_count, 0)
  expect_match(phyto$ProviderDiagnostics$message, "local standardized index")
})

test_that("resolver uses manifest-backed LOTUS lookup directories", {
  lotus = data.frame(
    traditional_name = c("Salicin", "Salixalbin"),
    lotus_id = c("LTS000010", "LTS000012"),
    allTaxa = c("Plantae; Salicaceae; Salix nigra",
                "Plantae; Salicaceae; Salix alba"),
    doi = c("10.1000/salix", "10.1000/alba"),
    stringsAsFactors = FALSE
  )
  index = standardizeLotusIndex(lotus)
  lookup_dir = tempfile()
  species_rows = cbind(
    index_key_type = "species", index_key = "salix nigra",
    index[index$species == "Salix nigra", , drop = FALSE]
  )
  genus_rows = cbind(
    index_key_type = "genus", index_key = "salix",
    index[index$species == "Salix alba", , drop = FALSE]
  )
  dir.create(file.path(lookup_dir, "keys", "species", "sa"),
             recursive = TRUE)
  dir.create(file.path(lookup_dir, "keys", "genus", "sa"),
             recursive = TRUE)
  utils::write.csv(species_rows,
                   file.path(lookup_dir, "keys", "species", "sa",
                             "salix_nigra.csv"),
                   row.names = FALSE, na = "")
  utils::write.csv(genus_rows,
                   file.path(lookup_dir, "keys", "genus", "sa", "salix.csv"),
                   row.names = FALSE, na = "")
  jsonlite::write_json(
    list(format = "uafR_lotus_lookup_index", version = 1,
         lookup_layout = "exact", shard_prefix_length = 2),
    file.path(lookup_dir, "manifest.json"),
    auto_unbox = TRUE
  )

  phyto = resolvePlantPhytochemistry(
    data.frame(species = "Salix nigra", family = "Salicaceae",
               stringsAsFactors = FALSE),
    sources = "lotus",
    lotus_index = lookup_dir,
    taxon_fallback = c("species", "genus"),
    enrich_compounds = FALSE,
    detail = "none"
  )
  occ = phyto$PlantCompoundOccurrences

  expect_equal(nrow(occ), 2)
  expect_true(any(occ$compound_name == "Salicin" &
                    occ$matched_rank == "species"))
  expect_true(any(occ$compound_name == "Salixalbin" &
                    occ$matched_rank == "genus"))
  expect_equal(phyto$ProviderDiagnostics$request_count, 0)
  expect_match(phyto$ProviderDiagnostics$message, "manifest-backed")
  expect_true(all(occ$analysis_ready[occ$matched_rank == "species"] == "Yes"))
  expect_true(all(occ$analysis_ready[occ$matched_rank == "genus"] == "No"))
})

test_that("buildLotusIndex writes a compact reusable LOTUS index", {
  lotus_a = data.frame(
    traditional_name = c("Salicin", "Salicin"),
    lotus_id = c("LTS000010", "LTS000010"),
    allTaxa = c("Plantae; Salicaceae; Salix nigra",
                "Plantae; Salicaceae; Salix nigra"),
    smiles = c("C1=CC=C(C=C1)CO", "C1=CC=C(C=C1)CO"),
    doi = c("10.1000/salix", "10.1000/salix"),
    stringsAsFactors = FALSE
  )
  lotus_b = data.frame(
    traditional_name = c("Salicin", "Achilleol A"),
    lotus_id = c("LTS000010", "LTS000011"),
    allTaxa = c("Plantae; Salicaceae; Salix nigra",
                "Plantae; Asteraceae; Achillea millefolium"),
    inchikey = c(NA_character_, "AAAAAAAAAAAAAA-UHFFFAOYSA-N"),
    pmid = c(NA_character_, "12345678"),
    doi = c("10.1000/salix", NA_character_),
    stringsAsFactors = FALSE
  )
  src_a = tempfile(fileext = ".csv")
  src_b = tempfile(fileext = ".csv")
  out_file = tempfile(fileext = ".csv")
  utils::write.csv(lotus_a, src_a, row.names = FALSE)
  utils::write.csv(lotus_b, src_b, row.names = FALSE)

  built = buildLotusIndex(
    input = c(src_a, src_b),
    out_file = out_file,
    overwrite = TRUE,
    progress = FALSE
  )
  manifest_file = sub("\\.csv$", "_manifest.json", out_file)
  reloaded = standardizeLotusIndex(out_file)
  occ = queryLotusIndex("Achillea millefolium", out_file)

  expect_s3_class(built, "uaf_lotus_index_build")
  expect_true(file.exists(out_file))
  expect_true(file.exists(manifest_file))
  expect_equal(nrow(built$BuildManifest), 2)
  expect_equal(built$BuildSummary$source_count, 2)
  expect_equal(built$BuildSummary$index_row_count, 2)
  expect_equal(built$BuildSummary$duplicate_row_count, 1)
  expect_equal(nrow(reloaded), 2)
  expect_equal(nrow(occ), 1)
  expect_equal(occ$compound_name, "Achilleol A")
})

test_that("buildLotusIndex records unreadable source diagnostics", {
  lotus = data.frame(
    traditional_name = "Salicin",
    allTaxa = "Plantae; Salicaceae; Salix nigra",
    stringsAsFactors = FALSE
  )
  src = tempfile(fileext = ".csv")
  bad = tempfile(fileext = ".zip")
  out_file = tempfile(fileext = ".rds")
  utils::write.csv(lotus, src, row.names = FALSE)
  writeLines("not a real zip", bad)

  built = buildLotusIndex(
    input = c(src, bad),
    out_file = out_file,
    format = "rds",
    overwrite = TRUE,
    progress = FALSE
  )

  expect_true(file.exists(out_file))
  expect_equal(nrow(built$LotusIndex), 1)
  expect_equal(sum(built$BuildManifest$status == "error"), 1)
  expect_match(
    built$BuildManifest$error_message[
      built$BuildManifest$status == "error"
    ],
    "zipped BSON dumps"
  )
})

test_that("standardizePlantCompoundIntake normalizes curated fallback rows", {
  intake = data.frame(
    species = "Salix nigra",
    compound_name = "Salicin",
    source_database = "manual",
    citation_or_url = "https://example.test/salix",
    evidence_tier = "manual_curated",
    plant_part = "bark",
    stringsAsFactors = FALSE
  )
  out = standardizePlantCompoundIntake(intake)
  expect_equal(nrow(out), 1)
  expect_equal(out$compound_name_clean, "salicin")
  expect_equal(out$confidence, "high")
  expect_equal(out$matched_rank, "species")
  expect_equal(out$occurrence_status, "curated_reported")
  expect_equal(out$plant_part_group, "bark_wood")
  expect_equal(out$analysis_ready, "Yes")
  expect_true(out$evidence_quality_score > 0.9)
})

test_that("filterPlantPhytochemistryEvidence keeps analysis-ready biological contexts", {
  intake = data.frame(
    species = c("Salix nigra", "Salix nigra", "Salix nigra"),
    compound_name = c("Salicin", "Caffeic acid", "Candidate"),
    source_database = c("manual", "manual", "PubTator"),
    citation_or_url = c("https://example.test/1",
                        "https://example.test/2",
                        "https://pubmed.ncbi.nlm.nih.gov/1/"),
    evidence_tier = c("manual_curated", "manual_curated",
                      "direct_species_pubtator_candidate"),
    plant_part = c("bark", "leaf", NA),
    stringsAsFactors = FALSE
  )
  occurrences = standardizePlantCompoundIntake(intake)
  bark = filterPlantPhytochemistryEvidence(
    occurrences,
    plant_part_group = "bark_wood"
  )
  all_evidence = filterPlantPhytochemistryEvidence(
    occurrences,
    occurrence_status = "all",
    analysis_ready = NULL,
    min_confidence = NULL
  )
  expect_equal(nrow(bark), 1)
  expect_equal(bark$compound_name_clean, "salicin")
  expect_true(all(bark$analysis_ready == "Yes"))
  expect_true(any(all_evidence$occurrence_status == "candidate"))
})

test_that("resolver collapses duplicate occurrence evidence keys", {
  curated = data.frame(
    species = c("Salix nigra", "Salix nigra"),
    compound_name = c("salicin", "salicin"),
    source_database = c("manual", "manual"),
    citation_or_url = c("https://example.test/salicin",
                        "https://example.test/salicin"),
    evidence_tier = c("manual_curated", "manual_curated"),
    plant_part = c("bark", NA),
    method = c(NA, "LC-MS"),
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = character(),
    curated_data = curated,
    enrich_compounds = FALSE,
    detail = "none"
  )
  expect_equal(nrow(phyto$PlantCompoundOccurrences), 1)
  expect_equal(phyto$PlantCompoundOccurrences$plant_part_group, "bark_wood")
  expect_equal(phyto$PlantCompoundOccurrences$method_group, "lc_ms")
  expect_false(any(phyto$Validation$Issues$issue ==
                     "Duplicate evidence keys were found"))
})

test_that("review table promotes candidate evidence only after review", {
  candidate = standardizePlantCompoundIntake(data.frame(
    species = "Salix nigra",
    compound_name = "salicin",
    source_database = "PubTator",
    citation_or_url = "https://pubmed.ncbi.nlm.nih.gov/333/",
    evidence_tier = "direct_species_pubtator_candidate",
    stringsAsFactors = FALSE
  ))
  review = plantPhytochemistryReviewTable(candidate)
  expect_equal(nrow(review), 1)
  expect_equal(review$review_decision, "needs_review")

  review$review_decision = "promote_curated"
  review$reviewed_by = "test reviewer"
  review$review_note = "Paper reports salicin in sampled bark."
  review$proposed_citation_or_url = "https://doi.org/10.1000/example"
  review$proposed_plant_part = "bark"
  review$proposed_method = "LC-MS"
  updated = applyPlantPhytochemistryReview(candidate, review)

  expect_equal(updated$occurrence_status, "curated_reported")
  expect_equal(updated$analysis_ready, "Yes")
  expect_equal(updated$plant_part_group, "bark_wood")
  expect_equal(updated$method_group, "lc_ms")
  expect_true(updated$evidence_quality_score > 0.9)

  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = character(),
    curated_data = candidate,
    enrich_compounds = FALSE,
    detail = "none"
  )
  phyto_review = plantPhytochemistryReviewTable(phyto)
  phyto_review$review_decision = "promote_curated"
  phyto_review$reviewed_by = "test reviewer"
  phyto_review$proposed_citation_or_url = "https://doi.org/10.1000/example"
  phyto_updated = applyPlantPhytochemistryReview(phyto, phyto_review)
  expect_s3_class(phyto_updated, "uaf_plant_phytochemistry")
  expect_equal(phyto_updated$PlantCompoundOccurrences$occurrence_status,
               "curated_reported")
  expect_true(any(phyto_updated$Provenance$source == "manual_review"))
})

test_that("resolver reports explicit no-hit diagnostics without fabricating compounds", {
  phyto = resolvePlantPhytochemistry("No hit plant",
                                     sources = "lotus",
                                     enrich_compounds = FALSE,
                                     detail = "none",
                                     request_fun = function(url) list(results = list()))
  expect_s3_class(phyto, "uaf_plant_phytochemistry")
  expect_equal(nrow(phyto$PlantCompoundOccurrences), 0)
  expect_equal(phyto$SpeciesChemistrySummary$query_status, "no_public_records")
  expect_true(any(phyto$ProviderDiagnostics$status == "no_records"))
})

test_that("provider diagnostics retain elapsed time and request errors", {
  phyto = resolvePlantPhytochemistry(
    "Salix nigra",
    sources = "lotus",
    enrich_compounds = FALSE,
    detail = "none",
    request_fun = function(url) stop("simulated provider timeout"),
    request_timeout = 1,
    cache = FALSE,
    throttle = 0,
    progress = FALSE
  )

  expect_equal(nrow(phyto$PlantCompoundOccurrences), 0)
  expect_equal(phyto$ProviderDiagnostics$status, "warning")
  expect_equal(phyto$ProviderDiagnostics$error_count, 1)
  expect_true("elapsed_seconds" %in% names(phyto$ProviderDiagnostics))
  expect_true("error_messages" %in% names(phyto$ProviderDiagnostics))
  expect_true(!is.na(phyto$ProviderDiagnostics$elapsed_seconds))
  expect_true(grepl("simulated provider timeout",
                    phyto$ProviderDiagnostics$error_messages))
})

test_that("provider text context does not mislabel methods as plant parts", {
  expect_true(is.na(.plant_provider_text_context_value(
    "gas chromatography and quadrupole mass spectrometry",
    "plant_part"
  )))
  expect_true(is.na(.plant_provider_text_context_value(
    "essential oil composition",
    "plant_part"
  )))
  expect_equal(.plant_provider_text_context_value(
    "leaf extract analyzed by GC-MS",
    "plant_part"
  ), "leaf")
  expect_equal(.plant_provider_text_context_value(
    "leaf extract analyzed by GC-MS",
    "method"
  ), "GC-MS")
  expect_equal(.plant_context_raw_method_like(c(
    "gas chromatography and quadrupole mass",
    "leaf"
  )), c(TRUE, FALSE))
  occurrence = standardizePlantCompoundIntake(data.frame(
    species = "Lavandula angustifolia",
    compound_name = "linalool",
    source_database = "manual",
    citation_or_url = "https://example.test",
    evidence_tier = "manual_curated",
    plant_part = "essential oil",
    method = "database record",
    stringsAsFactors = FALSE
  ))
  expect_true(is.na(occurrence$plant_part))
  expect_equal(occurrence$plant_part_group, "extract_unspecified")
  expect_true(is.na(occurrence$method))
  expect_equal(occurrence$method_group, "database_record")
  conflict = .plant_normalize_occurrences(data.frame(
    species = "Populus deltoides",
    compound_name = "quercetin",
    source_database = "PubTator",
    source_record_id = "123",
    evidence_tier = "direct_species_pubtator_candidate",
    plant_part = "gas chromatography and quadrupole mass",
    plant_part_group = "leaf",
    stringsAsFactors = FALSE
  ))
  expect_equal(conflict$plant_part, "leaf")
  expect_equal(conflict$plant_part_group, "leaf")
  context = data.frame(
    species = "Populus deltoides",
    compound_name = "quercetin",
    source_database = "PubTator",
    source_record_id = "123",
    context_type = "plant_part",
    raw_context_text = "leaves",
    normalized_context = "leaf",
    source_field = "evidence_text",
    extraction_rule = "pubtator_co_mention_sentence_regex:leaf_or_foliar",
    context_confidence = "low",
    evidence_basis = "provider_sentence_regex",
    stringsAsFactors = FALSE
  )
  applied = .plant_apply_context_evidence(conflict, context)
  expect_equal(applied$plant_part, "leaf")
  expect_equal(applied$plant_part_group, "leaf")
})

test_that("live provider adapters normalize mocked public responses", {
  provider_request = function(url) {
    if (grepl("lotus", url)) {
      return(list(results = list(list(
        name = "salicin",
        lotus_id = "LTS000001",
        organism = list(value = "Salix nigra",
                        cleaned_organism_id = "75706"),
        taxonomy = list(ncbi = list(cleaned_organism_id = "75706",
                                    organism_value = "Salix nigra",
                                    species = "Salix nigra")),
        reference = "PMID:12345678"
      ), list(
        name = "salicin",
        lotus_id = "LTS000002",
        organism = "Salix nigra",
        plant_part = "leaf",
        tissue = "leaf",
        method = "GC-MS",
        reference = "PMID:12345678"
      ))))
    }
    if (grepl("knapsack", url)) {
      return(paste(
        "<table>",
        "<tr><th>C_ID</th><th>Metabolite</th><th>Formula</th><th>Organism</th></tr>",
        paste0(
          "<tr><td><a href=\"information.php?word=C00000001\">",
          "C00000001</a></td><td>salicylic acid</td><td>C7H6O3</td>",
          "<td>Salix nigra</td></tr>"
        ),
        "</table>"
      ))
    }
    if (grepl("esearch.fcgi", url) && grepl("taxonomy", url)) {
      return(list(esearchresult = list(idlist = list("75706"))))
    }
    if (grepl("pug_view/data/taxonomy", url)) {
      return(list(Record = list(
        Reference = list(list(ReferenceNumber = 1,
                              SourceName = "PubChem",
                              URL = "https://pubchem.ncbi.nlm.nih.gov/taxonomy/75706")),
        Section = list(list(
          TOCHeading = "Chemicals and Bioactivities",
          Information = list(list(
            Name = "Associated compound",
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(
              String = "salicin",
              Markup = list(list(Start = 0, Length = 7,
                                 URL = "https://pubchem.ncbi.nlm.nih.gov/compound/439503"))
            )))
          )),
          Section = list(list(
            TOCHeading = "Metabolites",
            Description = "Metabolites reported for this organism.",
            Information = list(list(
              ReferenceNumber = 1,
              Value = list(ExternalTableName = paste0(
                "collection=consolidatedcompoundtaxonomy",
                "&view=concise_taxonomy",
                "&srccmpdkind=Metabolite"
              ))
            ))
          ))
        ))
      )))
    }
    if (grepl("sphinxql.cgi", url)) {
      return(list(SDQOutputSet = list(list(
        status = list(code = 0),
        totalCount = 1,
        rows = list(list(
          cid = "5146",
          cmpdname = "saligenin",
          synonym = "saligenin",
          taxname = "Salix nigra",
          srccmpdkind = "Metabolite",
          dsn = "KNApSAcK Species-Metabolite Database",
          pmids = "12345678",
          dois = "10.1000/test",
          evurls = "https://doi.org/10.1000/test",
          citations = paste(
            "Saligenin reported from Salix nigra leaf extract by LC-MS."
          )
        ))
      ))))
    }
    list()
  }

  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = c("lotus", "knapsack", "pubchem"),
    taxon_fallback = "species",
    enrich_compounds = FALSE,
    detail = "none",
    request_fun = provider_request,
    cache = FALSE,
    throttle = 0
  )

  expect_true(all(c("LOTUS", "KNApSAcK", "PubChem Taxonomy") %in%
                    phyto$PlantCompoundOccurrences$source_database))
  expect_true(all(phyto$ProviderDiagnostics$status == "ok"))
  expect_true(any(phyto$PlantCompoundOccurrences$compound_name == "salicylic acid"))
  expect_true(any(phyto$PlantCompoundOccurrences$compound_id == "439503"))
  expect_true(any(phyto$PlantCompoundOccurrences$compound_name == "saligenin"))
  expect_true(any(phyto$PlantCompoundOccurrences$occurrence_type ==
                    "pubchem_taxonomy_metabolite_table"))
  expect_true(any(phyto$PlantCompoundOccurrences$pmid == "12345678"))
  expect_true(any(phyto$PlantCompoundOccurrences$doi == "10.1000/test"))
  expect_true(all(phyto$PlantCompoundOccurrences$evidence_tier ==
                    "direct_species_database"))
  expect_true("elapsed_seconds" %in% names(phyto$ProviderDiagnostics))
  expect_true(all(!is.na(phyto$ProviderDiagnostics$elapsed_seconds)))
  lotus = phyto$PlantCompoundOccurrences[
    phyto$PlantCompoundOccurrences$source_database == "LOTUS", ,
    drop = FALSE
  ]
  expect_true(any(is.na(lotus$plant_part)))
  expect_true(any(lotus$plant_part == "leaf", na.rm = TRUE))
  expect_false(any(grepl("^[0-9]+$", lotus$plant_part)))
  knapsack = phyto$PlantCompoundOccurrences[
    phyto$PlantCompoundOccurrences$source_database == "KNApSAcK", ,
    drop = FALSE
  ]
  expect_true(any(grepl("information\\.php\\?word=C00000001",
                        knapsack$evidence_url)))
  expect_true(any(grepl("Metabolite=salicylic acid",
                        knapsack$evidence_text, fixed = TRUE)))
  expect_true(all(knapsack$biological_context_status == "context_missing"))
})

test_that("resolver combines mocked provider rows and conservative literature candidates", {
  lotus_rows = data.frame(
    species = "Salix nigra",
    family = "Salicaceae",
    compound_name = "salicin",
    source_database = "LOTUS",
    source_record_id = "LTS000001",
    evidence_text = "salicin reported in Salix nigra",
    evidence_url = "https://example.test/lotus/LTS000001",
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = c("lotus", "pubmed", "pubtator"),
    enrich_compounds = TRUE,
    detail = "research",
    provider_results = list(lotus = lotus_rows),
    enrichment_fun = plant_enrichment_fixture,
    request_fun = pubmed_fixture_request,
    pubtator_request_fun = pubtator_fixture_request,
    throttle = 0
  )

  expect_s3_class(phyto, "uaf_plant_phytochemistry")
  expect_true(any(phyto$PlantCompoundOccurrences$source_database == "LOTUS"))
  expect_true(any(phyto$PlantCompoundOccurrences$source_database == "PubTator"))
  expect_true(any(phyto$PlantCompoundOccurrences$evidence_tier ==
                    "direct_species_pubtator_candidate"))
  expect_false(any(phyto$PlantCompoundOccurrences$source_database == "PubTator" &
                     phyto$PlantCompoundOccurrences$evidence_tier ==
                       "direct_species_database"))
  expect_true(any(phyto$LiteratureCandidates$source_database == "PubMed"))
  expect_true(any(phyto$LiteratureCandidates$source_database == "PubTator"))
  expect_true(any(phyto$CompoundResolution$resolved))
  expect_true(phyto$SpeciesChemistrySummary$compound_count >= 1)
  expect_true("compound__salicin" %in% names(phyto$SpeciesChemistryMatrix))
  expect_equal(phyto$Validation$Summary$Status, "pass")
})

test_that("plant chemistry comparability separates analysis scopes", {
  curated = data.frame(
    species = c(rep("Salix nigra", 5), rep("Zea mays", 3)),
    compound_name = c("salicin", "caffeine", "limonene", "glucose",
                      "linoleic acid", "alliin", "abscisic acid",
                      "unknownoid"),
    source_database = "manual",
    citation_or_url = paste0("https://example.test/comparison/", 1:8),
    evidence_tier = "manual_curated",
    plant_part = c("bark", "leaf", "leaf", "leaf", "seed", "root",
                   "leaf", "leaf"),
    method = c("LC-MS", "LC-MS", "GC-MS", "LC-MS", "LC-MS", "LC-MS",
               "LC-MS", "LC-MS"),
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = c("Salix nigra", "Zea mays"),
    sources = character(),
    curated_data = curated,
    enrich_compounds = TRUE,
    detail = "research",
    enrichment_fun = plant_comparability_enrichment_fixture
  )
  comparability = plantChemistryComparability(phyto)

  salicin = comparability[comparability$compound_name_clean == "salicin", ]
  limonene = comparability[comparability$compound_name_clean == "limonene", ]
  glucose = comparability[comparability$compound_name_clean == "glucose", ]
  linoleic = comparability[
    comparability$compound_name_clean == "linoleic_acid", ]
  unknown = comparability[comparability$compound_name_clean == "unknownoid", ]

  expect_equal(salicin$comparison_scope, "specialized_metabolites")
  expect_equal(salicin$biosynthetic_family,
               "phenolic_phenylpropanoid")
  expect_equal(limonene$comparison_scope,
               "volatile_specialized_metabolites")
  expect_equal(limonene$chemical_behavior, "volatile_semivolatile")
  expect_equal(glucose$comparison_scope, "primary_metabolites")
  expect_equal(linoleic$comparison_scope, "lipids_fatty_acids")
  expect_equal(unknown$comparison_scope, "unknown")
  expect_equal(unknown$comparable_for_matrix, "No")

  specialized = plantComparableChemistryMatrix(
    phyto,
    comparison_scope = "specialized_metabolites",
    feature = "comparison_group",
    mode = "binary"
  )
  volatile = plantComparableChemistryMatrix(
    phyto,
    comparison_scope = "volatile_specialized_metabolites",
    feature = "comparison_group",
    mode = "binary"
  )
  primary = plantComparableChemistryMatrix(
    phyto,
    comparison_scope = "primary_metabolites",
    feature = "comparison_group",
    mode = "binary"
  )

  expect_true("comparison_group__specialized_metabolites__phenolic_phenylpropanoid" %in%
                names(specialized))
  expect_false(any(grepl("primary_metabolites", names(specialized),
                         fixed = TRUE)))
  expect_true("comparison_group__volatile_specialized_metabolites__volatile_terpenoid" %in%
                names(volatile))
  expect_true("comparison_group__primary_metabolites__carbohydrate" %in%
                names(primary))
  expect_true("ChemistryComparability" %in% names(phyto))
  expect_true("ComparableChemistryMatrix" %in% names(phyto))
  expect_true(phyto$SpeciesChemistrySummary$comparable_primary_compound_count[
    phyto$SpeciesChemistrySummary$species == "Salix nigra"
  ] > 0)
})

test_that("context evidence extracts plant part and method from source text", {
  curated = data.frame(
    species = c("Salix nigra", "Zea mays"),
    compound_name = c("limonene", "limonene"),
    source_database = "manual",
    citation_or_url = c("https://example.test/context/leaf",
                        "https://example.test/context/root"),
    evidence_tier = "manual_curated",
    evidence_note = c(
      "GC-MS analysis of leaf essential oil reported limonene.",
      "GC-MS analysis of root exudates reported limonene."
    ),
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = c("Salix nigra", "Zea mays"),
    sources = character(),
    curated_data = curated,
    enrich_compounds = FALSE,
    detail = "none"
  )
  context = phyto$PlantContextEvidence

  expect_true("PlantContextEvidence" %in% names(phyto))
  expect_true(all(c("leaf", "exudate_rhizosphere", "gc_ms") %in%
                    context$normalized_context))
  expect_false("literature_curation" %in% context$normalized_context)
  expect_equal(phyto$PlantCompoundOccurrences$plant_part_group,
               c("leaf", "exudate_rhizosphere"))
  expect_equal(phyto$PlantCompoundOccurrences$method_group,
               c("gc_ms", "gc_ms"))
  expect_equal(phyto$PlantCompoundOccurrences$biological_context_status,
               c("plant_part_and_method_known",
               "plant_part_and_method_known"))
})

test_that("source-backed context enrichment applies PMID/DOI text conservatively", {
  occurrence = .plant_normalize_occurrences(data.frame(
    species = "Salix nigra",
    genus = "Salix",
    compound_name = "salicin",
    source_database = "LOTUS",
    source_record_id = "LTS000010",
    doi = "10.1000/salix",
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  sources = data.frame(
    source_database = "PubMed",
    source_record_id = "12345",
    pmid = "12345",
    doi = "10.1000/salix",
    title = "Phytochemical constituents from Salix nigra leaves",
    abstract = paste(
      "Salicin was identified from leaves of Salix nigra.",
      "The extract was analyzed by LC-MS."
    ),
    evidence_url = "https://pubmed.ncbi.nlm.nih.gov/12345/",
    stringsAsFactors = FALSE
  )

  context = enrichPlantContextEvidence(
    occurrence,
    context_sources = sources,
    apply = FALSE
  )
  updated = enrichPlantContextEvidence(
    list(PlantCompoundOccurrences = occurrence,
         PlantContextEvidence = .uaf_empty_table(.plant_context_evidence_cols())),
    context_sources = sources
  )

  expect_true(any(context$normalized_context == "leaf"))
  expect_true(any(context$normalized_context == "lc_ms"))
  expect_true(any(grepl("source_context", context$extraction_rule)))
  expect_true(all(context$evidence_basis %in%
                    c("source_backed_species_compound_sentence",
                      "source_backed_species_chemical_context_sentence",
                      "source_backed_document_context_sentence")))
  expect_equal(updated$PlantCompoundOccurrences$plant_part_group, "leaf")
  expect_equal(updated$PlantCompoundOccurrences$method_group, "lc_ms")
})

test_that("source-backed context enrichment rejects unrelated source sentences", {
  occurrence = .plant_normalize_occurrences(data.frame(
    species = "Salix nigra",
    genus = "Salix",
    compound_name = "salicin",
    source_database = "LOTUS",
    source_record_id = "LTS000010",
    doi = "10.1000/mixed",
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  unrelated = data.frame(
    source_database = "PubMed",
    source_record_id = "12346",
    pmid = "12346",
    doi = "10.1000/mixed",
    title = "Chemical profiling of unrelated plants",
    abstract = "Zea mays leaf extract was analyzed by GC-MS.",
    evidence_url = "https://pubmed.ncbi.nlm.nih.gov/12346/",
    stringsAsFactors = FALSE
  )

  context = enrichPlantContextEvidence(
    occurrence,
    context_sources = unrelated,
    apply = FALSE
  )

  expect_equal(nrow(context), 0)
})

test_that("source-backed context enrichment can fetch mocked PubMed abstracts", {
  occurrence = .plant_normalize_occurrences(data.frame(
    species = "Zea mays",
    genus = "Zea",
    compound_name = "benzoxazolinone",
    source_database = "LOTUS",
    source_record_id = "LTS000020",
    doi = "10.2000/zea",
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  requested = character()
  request_fun = function(url) {
    requested <<- c(requested, url)
    if (grepl("esearch.fcgi", url, fixed = TRUE)) {
      return('{"esearchresult":{"idlist":["999999"]}}')
    }
    paste0(
      "<PubmedArticle><MedlineCitation><PMID>999999</PMID>",
      "<Article><ArticleTitle>Root metabolites of Zea mays</ArticleTitle>",
      "<Abstract><AbstractText>",
      "Benzoxazolinone was detected in Zea mays root exudates by GC-MS.",
      "</AbstractText></Abstract></Article></MedlineCitation>",
      "<PubmedData><ArticleIdList>",
      "<ArticleId IdType=\"doi\">10.2000/zea</ArticleId>",
      "</ArticleIdList></PubmedData></PubmedArticle>"
    )
  }

  context = enrichPlantContextEvidence(
    occurrence,
    fetch_pubmed = TRUE,
    cache = FALSE,
    request_fun = request_fun,
    max_sources = 1,
    apply = FALSE
  )

  expect_true(any(grepl("esearch.fcgi", requested, fixed = TRUE)))
  expect_true(any(grepl("efetch.fcgi", requested, fixed = TRUE)))
  expect_true(any(context$normalized_context == "exudate_rhizosphere"))
  expect_true(any(context$normalized_context == "gc_ms"))
  expect_equal(unique(context$pmid), "999999")
})

test_that("PubMed context candidates prioritize informative shared sources", {
  occurrences = .plant_normalize_occurrences(data.frame(
    species = c("Salix nigra", "Salix nigra", "Zea mays"),
    genus = c("Salix", "Salix", "Zea"),
    compound_name = c("salicin", "catechin", "benzoxazolinone"),
    source_database = "LOTUS",
    source_record_id = paste0("LTS", 1:3),
    doi = c("10.1000/shared", "10.1000/shared", "10.1000/narrow"),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  sources = data.frame(
    source_database = "PubMed",
    source_record_id = "777777",
    pmid = "777777",
    doi = "10.1000/context",
    title = "GC-MS phytochemical profile of Salix nigra leaves",
    abstract = "Salicin was detected in leaf extracts from Salix nigra.",
    evidence_tier = "direct_species_literature",
    confidence = "medium",
    stringsAsFactors = FALSE
  )

  candidates = .plant_pubmed_context_candidates(occurrences, sources)
  shared = candidates[candidates$key == "10.1000/shared", , drop = FALSE]
  narrow = candidates[candidates$key == "10.1000/narrow", , drop = FALSE]

  expect_equal(candidates$key[[1]], "777777")
  expect_gt(shared$priority_score, narrow$priority_score)
  expect_equal(shared$source_text_signal, "0")
  expect_equal(candidates$source_text_signal[[1]], "2")
})

test_that("PubMed context fetch resolves highest-priority DOI first", {
  occurrences = .plant_normalize_occurrences(data.frame(
    species = c("Salix nigra", "Zea mays", "Zea mays", "Zea mays"),
    genus = c("Salix", "Zea", "Zea", "Zea"),
    compound_name = c("salicin", "benzoxazolinone", "ferulic acid",
                      "p-coumaric acid"),
    source_database = "LOTUS",
    source_record_id = paste0("LTS", 10:13),
    doi = c("10.1000/low", rep("10.1000/high", 3)),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  esearch_urls = character()
  request_fun = function(url) {
    decoded = utils::URLdecode(url)
    if (grepl("esearch.fcgi", decoded, fixed = TRUE)) {
      esearch_urls <<- c(esearch_urls, decoded)
      if (grepl("10.1000/high", decoded, fixed = TRUE)) {
        return('{"esearchresult":{"idlist":["222222"]}}')
      }
      return('{"esearchresult":{"idlist":["111111"]}}')
    }
    paste0(
      "<PubmedArticle><MedlineCitation><PMID>222222</PMID>",
      "<Article><ArticleTitle>Root exudates of Zea mays</ArticleTitle>",
      "<Abstract><AbstractText>",
      "Benzoxazolinone was detected in Zea mays root exudates by GC-MS.",
      "</AbstractText></Abstract></Article></MedlineCitation>",
      "<PubmedData><ArticleIdList>",
      "<ArticleId IdType=\"doi\">10.1000/high</ArticleId>",
      "</ArticleIdList></PubmedData></PubmedArticle>"
    )
  }

  fetched = .plant_fetch_pubmed_context_sources(
    occurrences = occurrences,
    context_sources = NULL,
    cache = FALSE,
    cache_dir = tempdir(),
    throttle = 0,
    ncbi_email = "",
    ncbi_tool = "uafR",
    ncbi_api_key = "",
    max_sources = 1,
    request_fun = request_fun,
    request_timeout = 5
  )

  expect_match(esearch_urls[[1]], "10.1000/high")
  expect_equal(unique(fetched$pmid), "222222")
})

test_that("duplicate context evidence keys collapse deterministically", {
  duplicate_context = data.frame(
    species = "Salix nigra",
    genus = "Salix",
    family = NA_character_,
    compound_name = "salicin",
    compound_name_clean = "salicin",
    source_database = "PubMed",
    source_record_id = "12345",
    pmid = "12345",
    evidence_url = "https://pubmed.ncbi.nlm.nih.gov/12345/",
    context_type = "plant_part",
    raw_context_text = c("leaves", "leaf extract"),
    normalized_context = "leaf",
    source_field = "abstract",
    extraction_rule = c("source_context:pubmed:leaf_plural",
                        "source_context:pubmed:leaf_extract"),
    context_confidence = c("low", "medium"),
    evidence_basis = c("source_backed_document_context_sentence",
                       "source_backed_species_compound_sentence"),
    requires_review = c("Yes", "No"),
    retrieved_at = "2026-06-10T00:00:00-0400",
    stringsAsFactors = FALSE
  )

  normalized = .plant_normalize_context_evidence(duplicate_context)

  expect_equal(nrow(normalized), 1)
  expect_equal(normalized$context_confidence, "medium")
  expect_equal(normalized$requires_review, "No")
  expect_match(normalized$raw_context_text, "leaves")
  expect_match(normalized$raw_context_text, "leaf extract")
  expect_match(normalized$extraction_rule, "leaf_plural")
  expect_match(normalized$extraction_rule, "leaf_extract")
})

test_that("candidate context evidence is flagged for review", {
  candidate = data.frame(
    species = "Salix nigra",
    compound_name = "salicin",
    source_database = "PubTator",
    citation_or_url = "https://example.test/pubtator/1",
    evidence_tier = "direct_species_pubtator_candidate",
    evidence_note = paste(
      "A PubTator co-mention candidate links Salix nigra and salicin",
      "in a leaf extract study."
    ),
    stringsAsFactors = FALSE
  )
  occurrence = standardizePlantCompoundIntake(candidate)
  context = plantContextEvidence(occurrence)
  low_confidence = context$context_confidence == "low"

  expect_true(any(low_confidence))
  expect_true(all(context$requires_review[low_confidence] == "Yes"))
})

test_that("provider-specific context extraction is conservative and auditable", {
  provider_rows = data.frame(
    species = c("Salix nigra", "Salix nigra", "Salix nigra"),
    compound_name = c("salicin", "caffeic acid", "salicylic acid"),
    source_database = c("KNApSAcK", "PubTator", "PubTator"),
    source_record_id = c("C000001", "PMID1", "PMID2"),
    evidence_url = paste0("https://example.test/provider-context/", 1:3),
    evidence_text = c(
      "C000001 | salicin | Salix nigra bark extract | LC-MS",
      paste(
        "Salix nigra and caffeic acid were mentioned together.",
        "A separate Zea mays leaf extract was analyzed by GC-MS."
      ),
      "Salix nigra leaf extract contained salicylic acid by GC-MS."
    ),
    evidence_tier = c("direct_species_database",
                      "direct_species_pubtator_candidate",
                      "direct_species_pubtator_candidate"),
    confidence = c("high", "low", "low"),
    stringsAsFactors = FALSE
  )
  occurrences = .plant_normalize_occurrences(provider_rows)
  context = plantContextEvidence(occurrences)
  audit = plantProviderContextAudit(
    list(PlantCompoundOccurrences = occurrences,
         PlantContextEvidence = context)
  )
  knapsack = context[context$source_database == "KNApSAcK", , drop = FALSE]
  pubtator = context[context$source_database == "PubTator", , drop = FALSE]

  expect_true(any(knapsack$normalized_context == "bark_wood"))
  expect_true(any(knapsack$normalized_context == "lc_ms"))
  expect_true(any(grepl("knapsack_row_text_regex",
                        knapsack$extraction_rule)))
  expect_false(any(pubtator$compound_name == "caffeic acid" &
                     pubtator$normalized_context == "leaf"))
  expect_true(any(pubtator$compound_name == "salicylic acid" &
                    pubtator$normalized_context == "leaf"))
  expect_true(all(pubtator$requires_review == "Yes"))
  expect_true("ProviderContextAudit" %in%
                unique(plantPhytochemistrySchema()$Table))
  expect_true(all(c("source_database", "context_known_fraction",
                    "audit_status", "recommended_action") %in%
                    names(audit)))
  expect_true(any(audit$source_database == "KNApSAcK"))
  expect_true(any(audit$source_database == "PubTator"))
  expect_true(audit$context_evidence_count[
    audit$source_database == "PubTator"
  ] > 0)
})

test_that("provider context audit distinguishes sparse KNApSAcK provenance", {
  provider_rows = data.frame(
    species = "Salix nigra",
    compound_name = "salicylic acid",
    source_database = "KNApSAcK",
    source_record_id = "C00000001",
    evidence_url =
      "https://www.knapsackfamily.com/knapsack_core/information.php?word=C00000001",
    evidence_text = paste(
      "C_ID=C00000001; Metabolite=salicylic acid;",
      "Molecular_formula=C7H6O3; Organism=Salix nigra"
    ),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  occurrences = .plant_normalize_occurrences(provider_rows)
  audit = plantProviderContextAudit(
    list(PlantCompoundOccurrences = occurrences)
  )

  expect_equal(audit$audit_status, "context_sparse")
  expect_match(audit$recommended_action,
               "KNApSAcK organism rows usually support")
  expect_match(audit$recommended_action,
               "evidence_url/source records")
})

test_that("context-aware comparable matrices keep like biological context", {
  curated = data.frame(
    species = c("Salix nigra", "Zea mays"),
    compound_name = c("limonene", "limonene"),
    source_database = "manual",
    citation_or_url = c("https://example.test/context-matrix/leaf",
                        "https://example.test/context-matrix/root"),
    evidence_tier = "manual_curated",
    evidence_note = c(
      "GC-MS analysis of leaf essential oil reported limonene.",
      "GC-MS analysis of root exudates reported limonene."
    ),
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = c("Salix nigra", "Zea mays"),
    sources = character(),
    curated_data = curated,
    enrich_compounds = TRUE,
    detail = "research",
    enrichment_fun = plant_comparability_enrichment_fixture
  )
  leaf = plantComparableChemistryMatrix(
    phyto,
    comparison_scope = "volatile_specialized_metabolites",
    feature = "comparison_group",
    plant_part_group = "leaf",
    require_context = TRUE
  )
  exudate = plantComparableChemistryMatrix(
    phyto,
    comparison_scope = "volatile_specialized_metabolites",
    feature = "comparison_group",
    plant_part_group = "exudate_rhizosphere",
    require_context = TRUE
  )

  expect_equal(leaf$species, "Salix nigra")
  expect_equal(exudate$species, "Zea mays")
  expect_true(any(grepl("volatile_specialized_metabolites", names(leaf),
                         fixed = TRUE)))
  expect_true(any(grepl("volatile_specialized_metabolites", names(exudate),
                         fixed = TRUE)))
})

test_that("source-backed class fields drive plant comparability before names", {
  curated = data.frame(
    species = "Salix nigra",
    compound_name = c("plain phenolic", "plain terpene", "plain alkaloid",
                      "plain lipid", "plain polyketide", "glucose",
                      "plain volatile"),
    source_database = "manual",
    citation_or_url = paste0("https://example.test/rich-class/", 1:7),
    evidence_tier = "manual_curated",
    stringsAsFactors = FALSE
  )
  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = character(),
    curated_data = curated,
    enrich_compounds = TRUE,
    detail = "research",
    enrichment_fun = plant_rich_class_enrichment_fixture
  )
  comparability = plantChemistryComparability(phyto)
  row_for = function(compound) {
    comparability[comparability$compound_name_clean ==
                    .plant_clean_compound(compound), , drop = FALSE]
  }

  phenolic = row_for("plain phenolic")
  terpene = row_for("plain terpene")
  alkaloid = row_for("plain alkaloid")
  lipid = row_for("plain lipid")
  polyketide = row_for("plain polyketide")
  glucose = row_for("glucose")
  volatile = row_for("plain volatile")

  expect_equal(phenolic$biosynthetic_family,
               "phenolic_phenylpropanoid")
  expect_equal(phenolic$classification_source_table, "ChemicalClasses")
  expect_true(grepl("Shikimates and Phenylpropanoids",
                    phenolic$classification_source_value))
  expect_equal(phenolic$comparability_confidence, "high")

  expect_equal(terpene$biosynthetic_family, "terpenoid")
  expect_equal(terpene$classification_source_table,
               "PubChemClassifications")
  expect_true(grepl("Prenol lipids", terpene$classification_source_value))

  expect_equal(alkaloid$biosynthetic_family, "alkaloid_nitrogenous")
  expect_equal(alkaloid$classification_source_table, "LOTUSProfile")

  expect_equal(lipid$comparison_scope, "lipids_fatty_acids")
  expect_equal(lipid$classification_source_table, "KEGGClassifications")

  expect_equal(polyketide$biosynthetic_family, "polyketide")
  expect_equal(polyketide$classification_source_table,
               "ChemicalTraitOntology")

  expect_equal(glucose$biosynthetic_family, "alkaloid_nitrogenous")
  expect_equal(glucose$classification_source_table, "ChemicalClasses")
  expect_false(grepl("compound_name_rule", glucose$comparability_basis,
                     fixed = TRUE))

  expect_equal(volatile$comparison_scope,
               "volatile_specialized_metabolites")
  expect_equal(volatile$classification_source_table, "DerivedGroups")
  expect_equal(volatile$chemical_behavior, "volatile_semivolatile")
})

test_that("enrichPlantCompounds accepts mocked categorate-like output", {
  intake = data.frame(species = "Camellia sinensis",
                      compound_name = "caffeine",
                      source_database = "manual",
                      citation_or_url = "https://example.test",
                      evidence_tier = "manual_curated",
                      stringsAsFactors = FALSE)
  occurrence = standardizePlantCompoundIntake(intake)
  enriched = enrichPlantCompounds(occurrence,
                                  detail = "research",
                                  enrichment_fun = plant_enrichment_fixture)
  expect_equal(enriched$CompoundResolution$resolved, TRUE)
  expect_equal(nrow(enriched$TraitEvidence), 1)
})

test_that("PubChem-only enrichment fallback resolves compounds and plant trait matrices", {
  intake = data.frame(species = "Salix nigra",
                      compound_name = "salicin",
                      source_database = "manual",
                      citation_or_url = "https://example.test",
                      evidence_tier = "manual_curated",
                      stringsAsFactors = FALSE)
  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = character(),
    curated_data = intake,
    enrich_compounds = TRUE,
    detail = "research",
    pubchem_fun = plant_pubchem_profile_fixture
  )

  expect_identical(phyto$CategorateResult$EnrichmentMode, "pubchem_only")
  expect_true(all(phyto$CompoundResolution$resolved))
  expect_true(nrow(phyto$CategorateResult$ChemicalTraits) > 0)
  expect_true(any(grepl("^chemtrait__", names(phyto$SpeciesChemistryMatrix))))
  expect_equal(phyto$SpeciesChemistrySummary$resolved_compound_count, 1)
  expect_equal(phyto$SpeciesChemistrySummary$uafR_validation_status, "pass")
})

test_that("PubChem-only enrichment can run in resumable batches", {
  intake = data.frame(species = "Salix nigra",
                      compound_name = c("salicin", "caffeine"),
                      source_database = "manual",
                      citation_or_url = "https://example.test",
                      evidence_tier = "manual_curated",
                      stringsAsFactors = FALSE)
  cache_dir = tempfile("plant_batch_cache_")
  phyto = resolvePlantPhytochemistry(
    plants = "Salix nigra",
    sources = character(),
    curated_data = intake,
    enrich_compounds = TRUE,
    detail = "research",
    pubchem_fun = plant_pubchem_profile_fixture,
    cache = TRUE,
    cache_dir = cache_dir,
    enrichment_batch_size = 1,
    progress = FALSE
  )

  expect_identical(phyto$CategorateResult$EnrichmentMode,
                   "pubchem_only_batched")
  expect_equal(phyto$CategorateResult$BatchCount, 2)
  expect_true(all(phyto$CompoundResolution$resolved))
  expect_true(any(file.exists(list.files(
    file.path(cache_dir, "compound_enrichment",
              "pubchem_only_batches"),
    full.names = TRUE
  ))))
})

test_that("identity-only plant compound resolution is fast and schema-compatible", {
  intake = data.frame(species = "Salix nigra",
                      compound_name = c("salicin", "caffeine"),
                      source_database = "manual",
                      citation_or_url = "https://example.test",
                      evidence_tier = "manual_curated",
                      stringsAsFactors = FALSE)
  occurrence = standardizePlantCompoundIntake(intake)
  resolved = resolvePlantCompoundIdentities(
    occurrence,
    pubchem_fun = plant_pubchem_profile_fixture,
    batch_size = 1,
    cache = TRUE,
    cache_dir = tempfile("plant_identity_cache_"),
    progress = FALSE
  )

  expect_true(all(.plant_compound_resolution_cols() %in% names(resolved)))
  expect_true(all(resolved$resolved))
  expect_true(all(!is.na(resolved$CID)))
  expect_identical(attr(resolved, "CategorateResult")$EnrichmentMode,
                   "pubchem_identity_batched")
})

test_that("LOTUS source structures resolve before PubChem and flag bad labels", {
  intake = data.frame(
    species = c("Salix nigra", "Salix nigra"),
    compound_name = c("salicin", "chamomile"),
    source_database = "LOTUS",
    citation_or_url = "https://lotus.test",
    evidence_tier = "direct_species_database",
    source_record_id = c("LTS_SALICIN", "LTS_CHAMOMILE"),
    stringsAsFactors = FALSE
  )
  occurrence = standardizePlantCompoundIntake(intake)
  occurrence$compound_id = occurrence$source_record_id
  occurrence$compound_id_type = "LOTUS"
  lotus_index = data.frame(
    species = c("Salix nigra", "Salix nigra"),
    compound_name = c("salicin", "chamomile"),
    lotus_id = c("LTS_SALICIN", "LTS_CHAMOMILE"),
    smiles = c("SALICIN_SMILES", "CHAMOMILE_STRUCTURE_SMILES"),
    inchikey = c("SALICINKEY", "CHAMOMILEKEY"),
    molecular_formula = c("C13H18O7", "C15H10O5"),
    stringsAsFactors = FALSE
  )

  resolved = resolvePlantCompoundIdentities(
    occurrence,
    lotus_index = lotus_index,
    cache = FALSE,
    pubchem_fun = function(...) stop("PubChem should not be called")
  )

  expect_true(all(resolved$resolved))
  expect_true(all(resolved$resolution_source == "LOTUS_source_identity"))
  expect_equal(resolved$SMILES[
    resolved$compound_name_clean == "salicin"
  ], "SALICIN_SMILES")

  review = plantCompoundIdentityReviewTable(
    list(PlantCompoundOccurrences = occurrence,
         CompoundResolution = resolved)
  )
  expect_equal(nrow(review), 1)
  expect_equal(review$compound_name, "chamomile")
  expect_equal(review$identity_issue_type,
               "resolved_source_or_product_label")
})

test_that("source-ambiguous names can resolve by PubChem but remain auditable", {
  intake = data.frame(
    species = c("Salix nigra", "Salix nigra"),
    compound_name = c("quercetin", "quercetin"),
    source_database = "LOTUS",
    citation_or_url = "https://lotus.test",
    evidence_tier = "direct_species_database",
    source_record_id = c("LTS_QUERCETIN_A", "LTS_QUERCETIN_B"),
    stringsAsFactors = FALSE
  )
  occurrence = standardizePlantCompoundIntake(intake)
  occurrence$compound_id = occurrence$source_record_id
  occurrence$compound_id_type = "LOTUS"
  lotus_index = data.frame(
    species = c("Salix nigra", "Salix nigra"),
    compound_name = c("quercetin", "quercetin"),
    lotus_id = c("LTS_QUERCETIN_A", "LTS_QUERCETIN_B"),
    smiles = c("QUERCETIN_STRUCTURE_A", "QUERCETIN_STRUCTURE_B"),
    inchikey = c("QUERCETINKEYA", "QUERCETINKEYB"),
    molecular_formula = c("C15H10O7", "C15H10O7"),
    stringsAsFactors = FALSE
  )
  pubchem_identity_only = function(compounds, ...) {
    list(
      identity = data.frame(Query = compounds,
                            CID = 5280343,
                            MatchStatus = "resolved",
                            QueriedName = compounds,
                            SourceURL = "https://pubchem.test/quercetin",
                            stringsAsFactors = FALSE),
      properties = data.frame(),
      synonyms = data.frame()
    )
  }

  resolved = resolvePlantCompoundIdentities(
    occurrence,
    lotus_index = lotus_index,
    cache = FALSE,
    pubchem_fun = pubchem_identity_only
  )

  expect_equal(nrow(resolved), 1)
  expect_true(resolved$resolved)
  expect_equal(as.integer(resolved$CID), 5280343)
  expect_equal(resolved$resolution_source, "pubchemProfile_identity")
  expect_match(resolved$notes, "Multiple LOTUS source structures")

  review = plantCompoundIdentityReviewTable(
    list(PlantCompoundOccurrences = occurrence,
         CompoundResolution = resolved)
  )
  expect_equal(nrow(review), 1)
  expect_equal(review$identity_issue_type,
               "resolved_source_structure_ambiguous")
  expect_equal(review$recommended_decision, "review_record_level_structure")
})

test_that("compound identity review decisions can be reapplied reproducibly", {
  occurrence = standardizePlantCompoundIntake(data.frame(
    species = "Salix nigra",
    compound_name = c("ambiguousoid", "source oil"),
    source_database = "LOTUS",
    citation_or_url = "https://lotus.test/review",
    evidence_tier = "direct_species_database",
    stringsAsFactors = FALSE
  ))
  resolution = data.frame(
    compound_name = c("ambiguousoid", "source oil"),
    compound_name_clean = .plant_clean_compound(c("ambiguousoid",
                                                  "source oil")),
    query_count = c(1, 1),
    resolved = c(FALSE, TRUE),
    CID = c(NA_character_, NA_character_),
    InChIKey = c(NA_character_, "SOURCEKEY"),
    SMILES = c(NA_character_, "SOURCE_SMILES"),
    MolecularFormula = c(NA_character_, "C10H20"),
    resolution_source = c("source_identity_review_required",
                          "LOTUS_source_identity"),
    notes = c("Multiple LOTUS source structures map to this normalized key.",
              ""),
    stringsAsFactors = FALSE
  )
  review = plantCompoundIdentityReviewTable(
    list(PlantCompoundOccurrences = occurrence,
         CompoundResolution = resolution),
    include_resolved = TRUE
  )
  review$review_decision[review$compound_name == "ambiguousoid"] =
    "update_identity"
  review$proposed_compound_name[review$compound_name == "ambiguousoid"] =
    "reviewed ambiguousoid"
  review$proposed_smiles[review$compound_name == "ambiguousoid"] =
    "REVIEWED_SMILES"
  review$proposed_inchikey[review$compound_name == "ambiguousoid"] =
    "REVIEWEDKEY"
  review$proposed_molecular_formula[review$compound_name == "ambiguousoid"] =
    "C9H10O2"
  review$proposed_resolution_source[review$compound_name == "ambiguousoid"] =
    "manual_identity_review"
  review$review_decision[review$compound_name == "source oil"] = "reject"
  review$reviewed_by = "test reviewer"

  updated = applyPlantCompoundIdentityReview(resolution, review)

  accepted = updated[updated$compound_name_clean ==
                       .plant_clean_compound("ambiguousoid"), ]
  rejected = updated[updated$compound_name_clean ==
                       .plant_clean_compound("source oil"), ]
  expect_true(accepted$resolved)
  expect_equal(accepted$compound_name, "reviewed ambiguousoid")
  expect_equal(accepted$SMILES, "REVIEWED_SMILES")
  expect_equal(accepted$resolution_source, "manual_identity_review")
  expect_false(rejected$resolved)
  expect_equal(rejected$resolution_source, "review_excluded")
})

test_that("compound resolution preserves PubChem SMILES fallback fields", {
  occurrence = standardizePlantCompoundIntake(data.frame(
    species = "Salix nigra",
    compound_name = "1-octen-3-ol",
    source_database = "manual",
    citation_or_url = "https://example.test/smiles",
    evidence_tier = "manual_curated",
    stringsAsFactors = FALSE
  ))
  categorate_result = list(
    PubChemProperties = data.frame(
      Query = "1-octen-3-ol",
      CID = 18827,
      MolecularFormula = "C8H16O",
      InChIKey = "VSMOENVRRABVKN-UHFFFAOYSA-N",
      CanonicalSMILES = NA_character_,
      IsomericSMILES = NA_character_,
      SMILES = "CCCCCC(C=C)O",
      ConnectivitySMILES = "CCCCCC(C=C)O",
      stringsAsFactors = FALSE
    ),
    EnrichmentMode = "pubchem_identity"
  )

  resolved = .plant_compound_resolution(occurrence, categorate_result)

  expect_true(resolved$resolved)
  expect_equal(resolved$SMILES, "CCCCCC(C=C)O")
})

test_that("compound keys preserve stereochemistry and Greek-letter variants", {
  compounds = c(paste0(intToUtf8(0x03b1), "-carotene"),
                paste0(intToUtf8(0x03b2), "-carotene"),
                "(+)-limonene", "(-)-limonene",
                paste(intToUtf8(0x03b2), "alanine"), "alanine")

  expect_equal(
    .plant_clean_compound(compounds),
    c("alpha_carotene", "beta_carotene", "plus_limonene",
      "minus_limonene", "beta_alanine", "alanine")
  )
})

test_that("compound resolution collapses punctuation aliases but not isomers", {
  occurrence = standardizePlantCompoundIntake(data.frame(
    species = "Salix nigra",
    compound_name = c("abscisic acid", "abscisic acid,",
                      paste0(intToUtf8(0x03b1), "-carotene"),
                      paste0(intToUtf8(0x03b2), "-carotene")),
    source_database = "manual",
    citation_or_url = "https://example.test/key-collapse",
    evidence_tier = "manual_curated",
    stringsAsFactors = FALSE
  ))
  categorate_result = list(
    PubChemProperties = data.frame(
      Query = c("abscisic acid", "abscisic acid,",
                paste0(intToUtf8(0x03b1), "-carotene"),
                paste0(intToUtf8(0x03b2), "-carotene")),
      CID = c(643732, 643732, 4369188, 5280489),
      MolecularFormula = c("C15H20O4", "C15H20O4", "C40H56", "C40H56"),
      InChIKey = c("ABAKEY", "ABAKEY", "ALPHAKEY", "BETAKEY"),
      CanonicalSMILES = NA_character_,
      IsomericSMILES = NA_character_,
      SMILES = c("ABA", "ABA", "ALPHA", "BETA"),
      ConnectivitySMILES = c("ABA", "ABA", "ALPHA", "BETA"),
      stringsAsFactors = FALSE
    ),
    EnrichmentMode = "pubchem_identity"
  )

  resolved = .plant_compound_resolution(occurrence, categorate_result)

  expect_equal(nrow(resolved), 3)
  expect_true(all(c("abscisic_acid", "alpha_carotene", "beta_carotene") %in%
                    resolved$compound_name_clean))
  expect_equal(resolved$SMILES[
    resolved$compound_name_clean == "alpha_carotene"
  ], "ALPHA")
  expect_equal(resolved$SMILES[
    resolved$compound_name_clean == "beta_carotene"
  ], "BETA")
  expect_match(resolved$notes[
    resolved$compound_name_clean == "abscisic_acid"
  ], "Collapsed aliases")
})

test_that("plant phytochemistry pilot panel is a reusable benchmark input", {
  panel = plantPhytochemistryPilotPanel()
  remediation = plantPhytochemistryPilotPanel("remediation")

  expect_true(all(c("species", "panel_role", "expected_data_depth",
                    "review_focus", "rationale") %in% names(panel)))
  expect_equal(nrow(panel), 15)
  expect_true("Salix nigra" %in% panel$species)
  expect_true(all(remediation$species %in% panel$species))
  expect_true(all(nzchar(panel$review_focus)))
})

test_that("chunked plant batch runner writes production artifacts", {
  provider_rows = data.frame(
    species = c("Salix nigra", "Camellia sinensis", "Zea mays"),
    compound_name = c("salicin", "caffeine", "DIMBOA"),
    source_database = "LOTUS",
    source_record_id = c("LTS1", "LTS2", "LTS3"),
    evidence_url = paste0("https://example.test/", 1:3),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("plant_batch_out_")
  phyto = runPlantPhytochemistryBatch(
    plants = c("Salix nigra", "Camellia sinensis", "Zea mays"),
    sources = "lotus",
    provider_results = list(lotus = provider_rows),
    out_dir = out_dir,
    species_chunk_size = 1,
    compound_resolution_profile = "identity",
    pubchem_fun = plant_pubchem_profile_fixture,
    compound_batch_size = 2,
    cache = TRUE,
    throttle = 0,
    progress = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(phyto, "uaf_plant_phytochemistry")
  expect_equal(nrow(phyto$BatchChunkManifest), 3)
  expect_equal(phyto$BatchRunManifest$plant_count, 3)
  expect_equal(phyto$BatchRunManifest$compound_resolution_profile, "identity")
  expect_true(all(c("attempted_compound_count",
                    "attempted_unresolved_compound_count",
                    "not_attempted_compound_count") %in%
                    names(phyto$BatchRunManifest)))
  expect_true(all(phyto$CompoundResolution$resolved))
  expect_true(file.exists(file.path(out_dir, "all_occurrences.csv")))
  expect_true(file.exists(file.path(out_dir, "analysis_ready_occurrences.csv")))
  expect_true(file.exists(file.path(out_dir, "review_required_occurrences.csv")))
  expect_true(file.exists(file.path(out_dir, "provider_context_audit.csv")))
  expect_true(file.exists(file.path(out_dir, "compound_identity_resolution.csv")))
  expect_true(file.exists(file.path(out_dir, "run_manifest.json")))
  expect_true("BatchRunManifest" %in%
                exportPlantPhytochemistryWorkbook(
                  phyto,
                  tempfile("plant_batch_export_"),
                  format = "csv",
                  overwrite = TRUE
                )$Table)
})

test_that("plant phytochemistry pilot writes summaries and context matrices", {
  provider_rows = data.frame(
    species = c("Salix nigra", "Salix nigra", "Zea mays",
                "Camellia sinensis"),
    compound_name = c("limonene", "salicin", "glucose", "caffeine"),
    source_database = "LOTUS",
    source_record_id = c("LTS10", "LTS11", "LTS12", "LTS13"),
    evidence_url = paste0("https://example.test/pilot/", 1:4),
    evidence_text = c(
      "GC-MS analysis of leaf essential oil reported limonene.",
      "LC-MS analysis of bark extract reported salicin.",
      "LC-MS analysis of root exudates reported glucose.",
      "LC-MS analysis of leaf extract reported caffeine."
    ),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("plant_pilot_out_")
  phyto = runPlantPhytochemistryPilot(
    plants = c("Salix nigra", "Camellia sinensis", "Zea mays"),
    sources = "lotus",
    provider_results = list(lotus = provider_rows),
    out_dir = out_dir,
    species_chunk_size = 2,
    compound_resolution_profile = "research",
    enrichment_fun = plant_comparability_enrichment_fixture,
    cache = FALSE,
    throttle = 0,
    progress = FALSE,
    overwrite = TRUE
  )

  expect_s3_class(phyto, "uaf_plant_phytochemistry")
  expect_true("PilotSummary" %in% names(phyto))
  expect_true("PilotMatrices" %in% names(phyto))
  expect_true("PilotQAReport" %in% names(phyto))
  expect_true("ProviderContextAudit" %in% names(phyto))
  expect_true(all(c("pilot_status", "recommended_next_step",
                    "context_known_fraction",
                    "review_required_count") %in%
                    names(phyto$PilotSummary)))
  expect_true(all(c("check", "status", "recommendation") %in%
                    names(phyto$PilotQAReport)))
  expect_true(any(phyto$PilotSummary$pilot_status == "pilot_ready"))
  expect_true(file.exists(file.path(out_dir, "pilot_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "pilot_review_needed.csv")))
  expect_true(file.exists(file.path(out_dir, "pilot_qa_report.csv")))
  expect_true(file.exists(file.path(out_dir,
                                    "pilot_provider_context_audit.csv")))
  expect_true(file.exists(file.path(out_dir, "pilot_export_manifest.csv")))
  expect_true(file.exists(file.path(out_dir, "pilot_run_manifest.json")))
  expect_true(file.exists(file.path(
    out_dir, "matrix_specialized_metabolites.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir, "matrix_volatile_specialized_metabolites.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir, "matrix_leaf_associated_chemistry.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir, "matrix_root_exudate_associated_chemistry.csv"
  )))
  expect_true(file.exists(file.path(
    out_dir, "matrix_primary_metabolites.csv"
  )))
  expect_true("Zea mays" %in%
                phyto$PilotMatrices$matrix_root_exudate_associated_chemistry$species)
  expect_true("Salix nigra" %in%
                phyto$PilotMatrices$matrix_volatile_specialized_metabolites$species)
  expect_true("pilot_summary" %in% phyto$PilotExportManifest$artifact)
  expect_true("pilot_qa_report" %in% phyto$PilotExportManifest$artifact)
  expect_true("pilot_provider_context_audit" %in%
                phyto$PilotExportManifest$artifact)
  expect_true("matrix_primary_metabolites" %in%
                phyto$PilotExportManifest$artifact)
})

test_that("pilot QA distinguishes skipped identity resolution from failed resolution", {
  provider_rows = data.frame(
    species = c("Salix nigra", "Camellia sinensis"),
    compound_name = c("salicin", "caffeine"),
    source_database = "LOTUS",
    source_record_id = c("LTS20", "LTS21"),
    evidence_url = paste0("https://example.test/no-resolution/", 1:2),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  out_dir = tempfile("plant_pilot_no_resolution_")
  phyto = runPlantPhytochemistryPilot(
    plants = c("Salix nigra", "Camellia sinensis"),
    sources = "lotus",
    provider_results = list(lotus = provider_rows),
    out_dir = out_dir,
    species_chunk_size = 2,
    compound_resolution_profile = "none",
    cache = FALSE,
    throttle = 0,
    progress = FALSE
  )
  identity_qa = phyto$PilotQAReport[
    phyto$PilotQAReport$check == "identity_resolved_fraction", ,
    drop = FALSE
  ]

  expect_true(all(phyto$CompoundResolution$resolution_source ==
                    "not_attempted"))
  expect_true(all(phyto$CompoundResolution$notes ==
                    "Compound resolution was not requested."))
  expect_equal(phyto$BatchRunManifest$attempted_compound_count, 0)
  expect_equal(phyto$BatchRunManifest$not_attempted_compound_count, 2)
  expect_equal(identity_qa$status, "warning")
  expect_true(is.na(identity_qa$value))
  expect_match(identity_qa$details, "No compound identity resolution")
})

test_that("batch runner can cap identity resolution without dropping occurrences", {
  provider_rows = data.frame(
    species = c("Salix nigra", "Salix nigra", "Camellia sinensis"),
    compound_name = c("salicin", "caffeic acid", "caffeine"),
    source_database = "LOTUS",
    source_record_id = c("LTS1", "LTS2", "LTS3"),
    evidence_url = paste0("https://example.test/", 1:3),
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  )
  phyto = runPlantPhytochemistryBatch(
    plants = c("Salix nigra", "Camellia sinensis"),
    sources = "lotus",
    provider_results = list(lotus = provider_rows),
    species_chunk_size = 2,
    compound_resolution_profile = "identity",
    pubchem_fun = plant_pubchem_profile_fixture,
    max_unique_compounds = 1,
    cache = FALSE,
    throttle = 0,
    progress = FALSE
  )

  expect_equal(nrow(phyto$PlantCompoundOccurrences), 3)
  expect_equal(phyto$BatchRunManifest$resolution_input_occurrence_count, 1)
  expect_equal(phyto$BatchRunManifest$attempted_compound_count, 1)
  expect_equal(phyto$BatchRunManifest$not_attempted_compound_count, 2)
  expect_equal(sum(phyto$CompoundResolution$resolved), 1)
  expect_equal(sum(!(phyto$CompoundResolution$resolved %in% TRUE)), 2)
  expect_equal(sum(phyto$CompoundResolution$resolution_source ==
                     "not_attempted"), 2)
})

test_that("summary, matrix, joins, scores, and export work offline", {
  intake = data.frame(species = c("Zea mays", "Zea mays"),
                      compound_name = c("DIMBOA", "benzoxazolinone"),
                      source_database = "manual",
                      citation_or_url = "https://example.test",
                      evidence_tier = "manual_curated",
                      plant_part = c("root", "leaf"),
                      method = c("LC-MS", "GC-MS"),
                      stringsAsFactors = FALSE)
  phyto = resolvePlantPhytochemistry(
    plants = "Zea mays",
    sources = character(),
    curated_data = intake,
    enrich_compounds = TRUE,
    detail = "research",
    enrichment_fun = plant_enrichment_fixture
  )
  summary = summarizePlantPhytochemistry(phyto)
  matrix = plantPhytochemistryMatrix(phyto, mode = "count")
  joined = joinPlantChemistryMetadata(data.frame(species = "Zea mays",
                                                group = "crop"),
                                      phyto)
  scores = scorePlantChemistryCandidates(phyto)
  filtered = filterPlantPhytochemistryEvidence(phyto,
                                               plant_part_group = "root_belowground")
  tmp = tempfile("plant_phyto_export_")
  manifest = exportPlantPhytochemistryWorkbook(phyto, tmp,
                                               format = "csv",
                                               overwrite = TRUE)
  tmp_ready = tempfile("plant_phyto_export_ready_")
  ready_manifest = exportPlantPhytochemistryWorkbook(phyto, tmp_ready,
                                                     format = "csv",
                                                     preset = "analysis_ready",
                                                     overwrite = TRUE)
  manifest_file = manifest$FileName[manifest$Table == "ExportManifest"]

  expect_equal(summary$compound_count, 2)
  expect_equal(summary$analysis_ready_compound_count, 2)
  expect_true(summary$mean_evidence_quality_score > 0.9)
  expect_true(any(grepl("^compound__", names(matrix))))
  expect_true(any(grepl("^plant_part__", names(matrix))))
  expect_equal(nrow(filtered$PlantCompoundOccurrences), 1)
  expect_equal(filtered$PlantCompoundOccurrences$plant_part_group,
               "root_belowground")
  expect_equal(joined$group, "crop")
  expect_true(scores$chemistry_priority_score >= 0)
  expect_true(file.exists(file.path(tmp, manifest_file)))
  expect_true("PlantCompoundOccurrences" %in% manifest$Table)
  expect_true("PlantContextEvidence" %in% manifest$Table)
  expect_true("ProviderContextAudit" %in% manifest$Table)
  expect_true("ChemistryComparability" %in% manifest$Table)
  expect_true("ComparableChemistryMatrix" %in% manifest$Table)
  expect_true("PlantCompoundOccurrences" %in% ready_manifest$Table)
})

test_that("validation catches missing columns and duplicate evidence keys", {
  broken = list(
    PlantQueries = data.frame(query_id = "plant_1"),
    PlantNameResolution = .uaf_empty_table(.plant_name_resolution_cols()),
    ProviderDiagnostics = .uaf_empty_table(.plant_provider_diagnostic_cols()),
    PlantCompoundOccurrences = rbind(
      standardizePlantCompoundIntake(data.frame(
        species = "Salix nigra",
        compound_name = "salicin",
        source_database = "manual",
        citation_or_url = "https://example.test",
        evidence_tier = "manual_curated"
      )),
      standardizePlantCompoundIntake(data.frame(
        species = "Salix nigra",
        compound_name = "salicin",
        source_database = "manual",
        citation_or_url = "https://example.test",
        evidence_tier = "manual_curated"
      ))
	    ),
	    PlantContextEvidence = .uaf_empty_table(.plant_context_evidence_cols()),
	    ProviderContextAudit =
	      .uaf_empty_table(.plant_provider_context_audit_cols()),
	    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    CompoundResolution = .uaf_empty_table(.plant_compound_resolution_cols()),
    SpeciesChemistrySummary = .uaf_empty_table(.plant_summary_cols()),
    SpeciesChemistryMatrix = .uaf_empty_table(c("species")),
    ChemistryComparability = .uaf_empty_table(.plant_comparability_cols()),
    ComparableChemistryMatrix = .uaf_empty_table(c("species")),
    TraitEvidence = .uaf_empty_table(.plant_trait_evidence_cols()),
    Validation = list(Summary = data.frame(), TableQuality = data.frame(),
                      Issues = data.frame()),
    DataDictionary = plantPhytochemistrySchema(),
    Provenance = .uaf_empty_table(.plant_provenance_cols())
  )
  validation = validatePlantPhytochemistryResult(broken)
  expect_true(any(validation$Issues$issue == "Required column is missing"))
  expect_true(any(validation$Issues$issue == "Duplicate evidence keys were found"))
})
