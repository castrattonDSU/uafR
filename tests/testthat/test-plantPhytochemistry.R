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
                    "SpeciesChemistrySummary", "Provenance") %in%
                    schema$Table))
  occurrence = plantPhytochemistrySchema("PlantCompoundOccurrences")
  expect_true(all(c("species", "compound_name", "evidence_tier",
                    "confidence", "occurrence_status", "analysis_ready",
                    "plant_part_group", "evidence_quality_score") %in%
                    occurrence$Column))
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
        "<tr><td>C00000001</td><td>salicylic acid</td><td>C7H6O3</td><td>Salix nigra</td></tr>",
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
          ))
        ))
      )))
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
  expect_true(all(phyto$PlantCompoundOccurrences$evidence_tier ==
                    "direct_species_database"))
  lotus = phyto$PlantCompoundOccurrences[
    phyto$PlantCompoundOccurrences$source_database == "LOTUS", ,
    drop = FALSE
  ]
  expect_true(any(is.na(lotus$plant_part)))
  expect_true(any(lotus$plant_part == "leaf", na.rm = TRUE))
  expect_false(any(grepl("^[0-9]+$", lotus$plant_part)))
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
    LiteratureCandidates = .uaf_empty_table(.plant_literature_cols()),
    CompoundResolution = .uaf_empty_table(.plant_compound_resolution_cols()),
    SpeciesChemistrySummary = .uaf_empty_table(.plant_summary_cols()),
    SpeciesChemistryMatrix = .uaf_empty_table(c("species")),
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
