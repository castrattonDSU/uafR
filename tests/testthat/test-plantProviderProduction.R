write_npass_fixture = function(path, header, rows) {
  writeLines(c(paste(header, collapse = "\t"),
               vapply(rows, paste, collapse = "\t", FUN.VALUE = character(1))),
             path, useBytes = TRUE)
}

npass_fixture_files = function(root) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  paths = file.path(root, c("general.txt", "structure.txt", "pairs.txt",
                            "species.txt"))
  names(paths) = c("general", "structure", "pairs", "species")
  write_npass_fixture(
    paths[["general"]],
    c("np_id", "inchikey", "pref_name", "iupac_name", "pubchem_id",
      "molecular_formula"),
    list(
      c("NPC0001", "RYYVLZVUVIJVGH-UHFFFAOYSA-N", "caffeine",
        "1,3,7-trimethylpurine-2,6-dione", "2519", "C8H10N4O2"),
      c("NPC0002", "PFTAWBLQPZVEMU-UHFFFAOYSA-N", "epicatechin",
        "epicatechin", "72276", "C15H14O6"),
      c("NPC0003", "NGFMICBWJRZIBI-UHFFFAOYSA-N", "salicin", "salicin",
        "439503", "C13H18O7")
    )
  )
  write_npass_fixture(
    paths[["structure"]],
    c("np_id", "InChI", "InChIKey", "SMILES"),
    list(
      c("NPC0001", "InChI=1S/C8H10N4O2", "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
      c("NPC0002", "InChI=1S/C15H14O6", "PFTAWBLQPZVEMU-UHFFFAOYSA-N",
        "C1C(C(OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O"),
      c("NPC0003", "InChI=1S/C13H18O7", "NGFMICBWJRZIBI-UHFFFAOYSA-N",
        "C1=CC(=C(C=C1CO)O)OC2C(C(C(C(O2)CO)O)O)O")
    )
  )
  write_npass_fixture(
    paths[["species"]],
    c("org_id", "org_name", "org_tax_level", "org_tax_id", "species_name",
      "genus_name", "family_name"),
    list(
      c("ORG1", "Camellia sinensis", "species", "4442",
        "Camellia sinensis", "Camellia", "Theaceae"),
      c("ORG2", "Camellia japonica", "species", "4441",
        "Camellia japonica", "Camellia", "Theaceae"),
      c("ORG3", "Salix nigra", "species", "75706", "Salix nigra", "Salix",
        "Salicaceae")
    )
  )
  write_npass_fixture(
    paths[["pairs"]],
    c("src_org_record_id", "src_org_pair_id", "src_org_pair", "org_id",
      "np_id", "new_cp_found", "org_isolation_part",
      "org_collect_location", "org_collect_time", "ref_type", "ref_id",
      "ref_id_type", "ref_url"),
    list(
      c("REC1", "PAIR1", "Camellia sinensis--NPC0001", "ORG1", "NPC0001",
        "no", "leaf", "China", "2019", "article", "12345678", "PMID",
        "https://pubmed.ncbi.nlm.nih.gov/12345678/"),
      c("REC2", "PAIR2", "Camellia japonica--NPC0002", "ORG2", "NPC0002",
        "no", "flower", "Japan", "2020", "article", "10.1000/example",
        "DOI", "https://doi.org/10.1000/example"),
      c("REC3", "PAIR3", "Salix nigra--NPC0003", "ORG3", "NPC0003", "no",
        "bark", "USA", "2018", "article", "87654321", "PubMed",
        "https://pubmed.ncbi.nlm.nih.gov/87654321/")
    )
  )
  paths
}

test_that("NPASS index builds and preserves exact and fallback evidence", {
  root = tempfile("npass_fixture_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  files = npass_fixture_files(file.path(root, "raw"))
  index_dir = file.path(root, "index")
  build = buildNpassIndex(
    general_info = files[["general"]],
    structures = files[["structure"]],
    species_pairs = files[["pairs"]],
    species_info = files[["species"]],
    out_dir = index_dir,
    chunk_size = 2,
    progress = FALSE
  )

  expect_equal(build$BuildSummary$pair_row_count, 3)
  expect_equal(build$BuildSummary$index_version, "3")
  expect_equal(build$BuildSummary$usable_occurrence_row_count, 3)
  expect_true(file.exists(file.path(index_dir, "manifest.json")))
  expect_true(all(nzchar(build$SourceManifest$md5)))

  occurrences = queryNpassIndex(
    "Camellia sinensis", index_dir,
    taxon_fallback = c("species", "genus")
  )
  expect_equal(nrow(occurrences), 2)
  expect_equal(sort(unique(occurrences$matched_rank)), c("genus", "species"))
  direct = occurrences[occurrences$matched_rank == "species", , drop = FALSE]
  fallback = occurrences[occurrences$matched_rank == "genus", , drop = FALSE]
  expect_equal(direct$compound_name, "caffeine")
  expect_equal(direct$pmid, "12345678")
  expect_equal(direct$plant_part_group, "leaf")
  expect_equal(direct$occurrence_status, "direct_reported")
  expect_equal(fallback$compound_name, "epicatechin")
  expect_equal(fallback$occurrence_status, "taxon_fallback")

  identity = attr(occurrences, "SourceCompoundIdentity")
  expect_s3_class(identity, "data.frame")
  expect_equal(nrow(identity), 2)
  expect_true(all(nzchar(identity$InChIKey)))
  expect_true(all(nzchar(identity$SMILES)))
  expect_equal(sort(identity$MolecularFormula), c("C15H14O6", "C8H10N4O2"))
})

test_that("NPASS resolver reports accounting, resources, and source identity", {
  root = tempfile("npass_resolver_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  files = npass_fixture_files(file.path(root, "raw"))
  index_dir = file.path(root, "index")
  buildNpassIndex(files[["general"]], files[["structure"]], files[["pairs"]],
                  files[["species"]], index_dir, chunk_size = 1,
                  progress = FALSE)

  result = resolvePlantPhytochemistry(
    plants = c("Camellia sinensis", "Unknown plant"),
    sources = "npass",
    taxon_fallback = "species",
    enrich_compounds = FALSE,
    detail = "none",
    provider_indexes = list(npass = index_dir),
    require_all_providers = TRUE,
    progress = FALSE
  )

  expect_equal(nrow(result$ProviderQueryAccounting), 2)
  expect_equal(result$ProviderQueryAccounting$query_status,
               c("records", "no_records"))
  expect_equal(result$ProviderResourceManifest$availability_status,
               "available")
  expect_match(result$ProviderResourceManifest$md5, "^[a-f0-9]{32}$")
  expect_match(result$ProviderResourceManifest$sha256, "^[a-f0-9]{64}$")
  expect_equal(nrow(result$SourceCompoundIdentity), 1)
  expect_equal(result$SourceCompoundIdentity$CID, 2519)
  expect_equal(result$Validation$Summary$ErrorCount, 0)
  expect_equal(result$Validation$Summary$WarningCount, 1)
})

test_that("provider resource SHA-256 hashing is platform independent", {
  path = tempfile("uafr_sha256_")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  writeBin(charToRaw("abc"), path)

  expect_identical(
    .plant_sha256_file(path),
    paste0(
      "ba7816bf8f01cfea414140de5dae2223",
      "b00361a396177a9cb410ff61f20015ad"
    )
  )
})

test_that("provider merge keeps one accounting row per plant and provider", {
  plant = "Salix nigra"
  lotus = resolvePlantPhytochemistry(
    plant, sources = "lotus", enrich_compounds = FALSE, detail = "none",
    provider_results = list(lotus = data.frame(
      species = plant, compound_name = "salicin", source_database = "LOTUS",
      source_record_id = "LTS1", evidence_tier = "direct_species_database",
      confidence = "high", stringsAsFactors = FALSE
    )), progress = FALSE
  )
  pubmed = resolvePlantPhytochemistry(
    plant, sources = "pubmed", enrich_compounds = FALSE, detail = "none",
    provider_results = list(pubmed = list(LiteratureCandidates = data.frame(
      species = plant, source_database = "PubMed", source_record_id = "123",
      pmid = "123", title = "Salicin profiling", evidence_tier =
        "direct_species_literature", confidence = "medium",
      stringsAsFactors = FALSE
    ))), progress = FALSE
  )
  merged = mergePlantPhytochemistryResults(lotus, pubmed)
  expect_equal(nrow(merged$ProviderQueryAccounting), 2)
  expect_equal(sort(merged$ProviderQueryAccounting$provider),
               c("lotus", "pubmed"))
  expect_equal(nrow(merged$PlantCompoundOccurrences), 1)
  expect_equal(nrow(merged$LiteratureCandidates), 1)
  expect_equal(merged$Validation$Summary$ErrorCount, 0)
})

test_that("PubMed readiness probes a valid EInfo utility endpoint", {
  requested = character()
  availability = plantProviderAvailability(
    sources = "pubmed", probe_live = TRUE,
    request_fun = function(url) {
      requested <<- c(requested, url)
      "{}"
    }
  )

  expect_equal(availability$availability_status, "available")
  expect_match(requested, "einfo[.]fcgi[?]db=pubmed&retmode=json$")
  expect_equal(
    availability$source_url,
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  )
})

test_that("KNApSAcK follows documented redirect and rejects infraspecific rows", {
  requests = character()
  request = function(url) {
    requests <<- c(requests, url)
    if (grepl("info[.]php", url)) {
      return("<script>location.href='result.php?sname=organism&word=Camellia%20sinensis'</script>")
    }
    paste(
      "<table>",
      "<tr><th>C_ID</th><th>CAS ID</th><th>Metabolite</th><th>Formula</th><th>Mw</th><th>Organism</th></tr>",
      "<tr><td>C00000001</td><td>58-08-2</td><td>caffeine</td><td>C8H10N4O2</td><td>194.19</td><td>Camellia sinensis (L.) Kuntze</td></tr>",
      "<tr><td>C00000002</td><td></td><td>variant compound</td><td>C1H2</td><td>14</td><td>Camellia sinensis var. assamica</td></tr>",
      "<tr><td>C00000003</td><td></td><td>wrong species</td><td>C2H4</td><td>28</td><td>Camellia japonica</td></tr>",
      "</table>"
    )
  }
  result = resolvePlantPhytochemistry(
    "Camellia sinensis", sources = "knapsack", taxon_fallback = "species",
    enrich_compounds = FALSE, detail = "none", request_fun = request,
    cache = FALSE, throttle = 0, progress = FALSE
  )
  expect_equal(nrow(result$PlantCompoundOccurrences), 1)
  expect_equal(result$PlantCompoundOccurrences$compound_name, "caffeine")
  expect_equal(result$PlantCompoundOccurrences$matched_taxon,
               "Camellia sinensis (L.) Kuntze")
  expect_true(any(grepl("info[.]php", requests)))
  expect_true(any(grepl("result[.]php", requests)))
})

test_that("verified aliases authorize exact provider matches", {
  plants = data.frame(
    original_species = "Camellia typo",
    species = "Camellia sinensis",
    verified_synonym = "Thea sinensis",
    stringsAsFactors = FALSE
  )
  lotus_index = data.frame(
    species = "Thea sinensis",
    genus = "Thea",
    family = "Theaceae",
    compound_name = "caffeine",
    lotus_id = "LTS_ALIAS_1",
    cid = "2519",
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    molecular_formula = "C8H10N4O2",
    stringsAsFactors = FALSE
  )
  lotus = queryLotusIndex(plants, lotus_index, taxon_fallback = "species")
  expect_equal(nrow(lotus), 1)
  expect_equal(lotus$species, "Camellia sinensis")
  expect_equal(lotus$matched_taxon, "Thea sinensis")
  expect_equal(lotus$matched_rank, "species")
  expect_equal(lotus$evidence_tier, "direct_species_database")
  identity = attr(lotus, "SourceCompoundIdentity")
  expect_equal(identity$CID, 2519)
  expect_equal(identity$InChIKey, "RYYVLZVUVIJVGH-UHFFFAOYSA-N")

  requested = character()
  request = function(url) {
    requested <<- c(requested, utils::URLdecode(url))
    if (grepl("Thea sinensis", utils::URLdecode(url), fixed = TRUE)) {
      return(paste(
        "<table>",
        "<tr><th>C_ID</th><th>CAS ID</th><th>Metabolite</th><th>Formula</th><th>Mw</th><th>Organism</th></tr>",
        "<tr><td>C00000001</td><td>58-08-2</td><td>caffeine</td><td>C8H10N4O2</td><td>194.19</td><td>Thea sinensis (L.) Kuntze</td></tr>",
        "</table>"
      ))
    }
    "<table></table>"
  }
  knapsack = resolvePlantPhytochemistry(
    plants, sources = "knapsack", taxon_fallback = "species",
    enrich_compounds = FALSE, detail = "none", request_fun = request,
    cache = FALSE, throttle = 0, progress = FALSE
  )
  expect_equal(nrow(knapsack$PlantCompoundOccurrences), 1)
  expect_equal(knapsack$PlantCompoundOccurrences$species,
               "Camellia sinensis")
  expect_equal(knapsack$PlantCompoundOccurrences$matched_taxon,
               "Thea sinensis (L.) Kuntze")
  expect_true(any(grepl("Thea sinensis", requested, fixed = TRUE)))
  expect_false(any(grepl("Camellia typo", requested, fixed = TRUE)))
  expect_equal(knapsack$SourceCompoundIdentity$source_compound_id,
               "C00000001")
  expect_equal(knapsack$SourceCompoundIdentity$MolecularFormula,
               "C8H10N4O2")
})

test_that("PubChem taxonomy rejects mismatched TaxIDs and accepts verified aliases", {
  plants = data.frame(
    species = "Camellia sinensis",
    verified_synonym = "Thea sinensis",
    stringsAsFactors = FALSE
  )
  requested = character()
  request = function(url) {
    decoded = utils::URLdecode(url)
    requested <<- c(requested, decoded)
    if (grepl("esearch.fcgi", decoded, fixed = TRUE) &&
        grepl("Camellia sinensis", decoded, fixed = TRUE)) {
      return(list(esearchresult = list(idlist = list("999"))))
    }
    if (grepl("esearch.fcgi", decoded, fixed = TRUE) &&
        grepl("Thea sinensis", decoded, fixed = TRUE)) {
      return(list(esearchresult = list(idlist = list("4442"))))
    }
    if (grepl("esummary.fcgi", decoded, fixed = TRUE) &&
        grepl("id=999", decoded, fixed = TRUE)) {
      return(list(result = list(
        `999` = list(scientificname = "Camellia japonica", rank = "species")
      )))
    }
    if (grepl("esummary.fcgi", decoded, fixed = TRUE) &&
        grepl("id=4442", decoded, fixed = TRUE)) {
      return(list(result = list(
        `4442` = list(scientificname = "Thea sinensis", rank = "species")
      )))
    }
    if (grepl("pug_view/data/taxonomy/4442", decoded, fixed = TRUE)) {
      return(list(Record = list(
        Reference = list(list(
          ReferenceNumber = 1, SourceName = "PubChem",
          URL = "https://pubchem.ncbi.nlm.nih.gov/taxonomy/4442"
        )),
        Section = list(list(
          TOCHeading = "Chemicals and Bioactivities",
          Information = list(list(
            Name = "Associated compound", ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(list(
              String = "caffeine",
              Markup = list(list(
                Start = 0, Length = 8,
                URL = "https://pubchem.ncbi.nlm.nih.gov/compound/2519"
              ))
            )))
          ))
        ))
      )))
    }
    stop("Unexpected taxonomy fixture URL: ", decoded)
  }
  result = resolvePlantPhytochemistry(
    plants, sources = "pubchem", taxon_fallback = "species",
    enrich_compounds = FALSE, detail = "none", request_fun = request,
    cache = FALSE, throttle = 0, progress = FALSE
  )
  expect_equal(nrow(result$PlantCompoundOccurrences), 1)
  expect_equal(result$PlantCompoundOccurrences$matched_taxon,
               "Thea sinensis")
  expect_equal(result$PlantCompoundOccurrences$compound_id, "2519")
  expect_false(any(grepl("pug_view/data/taxonomy/999", requested,
                         fixed = TRUE)))
  expect_true(any(grepl("pug_view/data/taxonomy/4442", requested,
                        fixed = TRUE)))
  expect_equal(result$SourceCompoundIdentity$CID, 2519)
})

test_that("PubMed abstracts and PubTator BioC annotations remain linked", {
  pubmed_request = function(url) {
    if (grepl("esearch[.]fcgi", url)) {
      return(list(esearchresult = list(count = "7",
                                       idlist = list("111", "222"))))
    }
    if (grepl("esummary[.]fcgi", url)) {
      return(list(result = list(
        `111` = list(title = "Tea leaf chemical profile", articleids = list()),
        `222` = list(title = "Tea phytochemistry review", articleids = list())
      )))
    }
    if (grepl("efetch[.]fcgi", url)) {
      return(paste0(
        "<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>111</PMID>",
        "<Article><Abstract><AbstractText Label=\"RESULTS\">Caffeine was ",
        "measured in Camellia sinensis leaves by HPLC.</AbstractText></Abstract>",
        "</Article></MedlineCitation></PubmedArticle>",
        "<PubmedArticle><MedlineCitation><PMID>222</PMID><Article><Abstract>",
        "<AbstractText>Camellia sinensis phytochemistry was reviewed.</AbstractText>",
        "</Abstract></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
      ))
    }
    stop("Unexpected PubMed fixture URL: ", url)
  }
  pubtator_request = function(url) {
    expect_match(url, "publications/export/biocjson", fixed = TRUE)
    expect_match(url, "pmids=111,222", fixed = TRUE)
    expect_false(grepl("%2C", url, ignore.case = TRUE))
    list(documents = list(list(
      id = "111",
      passages = list(
        list(infons = list(type = "title"),
             text = "Tea leaf chemical profile", annotations = list()),
        list(infons = list(type = "abstract"),
             text = paste("Caffeine was measured in Camellia sinensis leaves",
                          "by HPLC."),
             annotations = list(
               list(text = "Caffeine", infons = list(type = "Chemical")),
               list(text = "Camellia sinensis",
                    infons = list(type = "Species"))
             ))
      )
    )))
  }

  result = resolvePlantPhytochemistry(
    "Camellia sinensis", sources = c("pubmed", "pubtator"),
    enrich_compounds = FALSE, detail = "none", max_pubmed_records = 2,
    request_fun = pubmed_request, pubtator_request_fun = pubtator_request,
    cache = FALSE, throttle = 0, progress = FALSE
  )

  pubmed = result$LiteratureCandidates[
    result$LiteratureCandidates$source_database == "PubMed", , drop = FALSE
  ]
  pubtator = result$LiteratureCandidates[
    result$LiteratureCandidates$source_database == "PubTator", , drop = FALSE
  ]
  expect_equal(nrow(pubmed), 2)
  expect_true(all(pubmed$literature_total_hit_count == 7))
  expect_true(all(pubmed$literature_search_truncated == "Yes"))
  expect_true(all(pubmed$abstract_retrieval_status == "retrieved"))
  expect_match(pubmed$abstract[pubmed$pmid == "111"], "measured")
  expect_equal(nrow(pubtator), 1)
  expect_equal(pubtator$chemical_mention, "Caffeine")
  expect_match(pubtator$evidence_text, "Camellia sinensis")
  expect_true(any(result$PlantCompoundOccurrences$occurrence_status ==
                    "candidate"))
  accounting = result$ProviderQueryAccounting[
    result$ProviderQueryAccounting$provider == "pubmed", , drop = FALSE
  ]
  expect_equal(accounting$total_hit_count, 7)
  expect_equal(accounting$query_truncated, "Yes")
  expect_equal(result$ProviderDiagnostics$request_count, c(3L, 1L))
})

test_that("unindexed PubTator PMID batches are cached no-record results", {
  root = tempfile("pubtator_no_record_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root)
  queries = uafR:::.plant_queries("Camellia sinensis", "species")
  literature = data.frame(
    species = "Camellia sinensis", pmid = c("42458239", "42452147"),
    source_database = "PubMed", stringsAsFactors = FALSE
  )
  requested = 0L
  first = uafR:::.plant_query_pubtator_literature(
    queries, cache = TRUE, cache_dir = root, throttle = 0,
    request_fun = function(url) {
      requested <<- requested + 1L
      stop("HTTP status was '400 Bad Request'")
    }, request_timeout = 30, pubmed_literature = literature
  )

  expect_equal(requested, 1L)
  expect_equal(first$ProviderDiagnostics$status, "no_records")
  expect_equal(first$ProviderDiagnostics$error_count, 0L)
  expect_match(first$ProviderDiagnostics$message,
               "2 PMID[(]s[)] had no indexed PubTator record")

  second = uafR:::.plant_query_pubtator_literature(
    queries, cache = TRUE, cache_dir = root, throttle = 0,
    request_fun = function(url) stop("cached no-record result was not reused"),
    request_timeout = 30, pubmed_literature = literature
  )
  expect_equal(second$ProviderDiagnostics$status, "no_records")
  expect_equal(second$ProviderDiagnostics$cache_hit_count, 1L)
  expect_equal(second$ProviderDiagnostics$error_count, 0L)
})

test_that("PubMed and PubTator retain synonym matches as candidate evidence", {
  plants = data.frame(
    species = "Camellia sinensis",
    verified_synonym = "Thea sinensis",
    stringsAsFactors = FALSE
  )
  searched = NA_character_
  pubmed_request = function(url) {
    decoded = utils::URLdecode(url)
    if (grepl("esearch[.]fcgi", decoded)) {
      searched <<- decoded
      return(list(esearchresult = list(count = "1", idlist = list("333"))))
    }
    if (grepl("esummary[.]fcgi", decoded)) {
      return(list(result = list(
        `333` = list(title = "Tea chemistry", articleids = list())
      )))
    }
    if (grepl("efetch[.]fcgi", decoded)) {
      return(paste0(
        "<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>333</PMID>",
        "<Article><Abstract><AbstractText>Caffeine was measured in Thea ",
        "sinensis leaves.</AbstractText></Abstract></Article></MedlineCitation>",
        "</PubmedArticle></PubmedArticleSet>"
      ))
    }
    stop("Unexpected PubMed URL: ", decoded)
  }
  pubtator_request = function(url) {
    list(documents = list(list(
      id = "333",
      passages = list(list(
        infons = list(type = "abstract"),
        text = "Caffeine was measured in Thea sinensis leaves.",
        annotations = list(
          list(text = "Caffeine", infons = list(type = "Chemical")),
          list(text = "Thea sinensis", infons = list(type = "Species"))
        )
      ))
    )))
  }
  result = resolvePlantPhytochemistry(
    plants, sources = c("pubmed", "pubtator"),
    enrich_compounds = FALSE, detail = "none", cache = FALSE,
    throttle = 0, request_fun = pubmed_request,
    pubtator_request_fun = pubtator_request, progress = FALSE
  )
  expect_match(searched, "Camellia sinensis", fixed = TRUE)
  expect_match(searched, "Thea sinensis", fixed = TRUE)
  candidate = result$PlantCompoundOccurrences[
    result$PlantCompoundOccurrences$source_database == "PubTator", ,
    drop = FALSE
  ]
  expect_equal(nrow(candidate), 1)
  expect_equal(candidate$species, "Camellia sinensis")
  expect_equal(candidate$compound_name, "Caffeine")
  expect_equal(candidate$occurrence_status, "candidate")
  expect_equal(candidate$analysis_ready, "No")
  pubtator = result$LiteratureCandidates[
    result$LiteratureCandidates$source_database == "PubTator", , drop = FALSE
  ]
  expect_match(pubtator$species_mention, "Thea sinensis", fixed = TRUE)
  expect_equal(pubtator$evidence_tier,
               "direct_species_pubtator_candidate")
})

test_that("provider diagnostics report cache hits and service failures", {
  cache_dir = tempfile("provider_cache_")
  on.exit(unlink(cache_dir, recursive = TRUE, force = TRUE), add = TRUE)
  request_count = 0L
  request = function(url) {
    request_count <<- request_count + 1L
    if (grepl("esearch[.]fcgi", url)) {
      return('{"esearchresult":{"count":"1","idlist":["111"]}}')
    }
    if (grepl("esummary[.]fcgi", url)) {
      return('{"result":{"111":{"title":"Plant chemistry","articleids":[]}}}')
    }
    if (grepl("efetch[.]fcgi", url)) {
      return(paste0(
        "<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>111</PMID>",
        "<Article><Abstract><AbstractText>Plant chemistry.</AbstractText>",
        "</Abstract></Article></MedlineCitation></PubmedArticle>",
        "</PubmedArticleSet>"
      ))
    }
    stop("Unexpected URL")
  }
  first = resolvePlantPhytochemistry(
    "Salix nigra", sources = "pubmed", enrich_compounds = FALSE,
    detail = "none", cache = TRUE, cache_dir = cache_dir,
    request_fun = request, throttle = 0, progress = FALSE
  )
  expect_equal(request_count, 3L)
  expect_equal(first$ProviderDiagnostics$request_count, 3L)
  expect_equal(first$ProviderDiagnostics$cache_hit_count, 0L)

  second = resolvePlantPhytochemistry(
    "Salix nigra", sources = "pubmed", enrich_compounds = FALSE,
    detail = "none", cache = TRUE, cache_dir = cache_dir,
    request_fun = function(url) stop("Cache-only rerun made a request"),
    throttle = 0, progress = FALSE
  )
  expect_equal(second$ProviderDiagnostics$request_count, 3L)
  expect_equal(second$ProviderDiagnostics$cache_hit_count, 3L)
  expect_equal(second$ProviderDiagnostics$status, "ok")

  busy = resolvePlantPhytochemistry(
    "Salix nigra", sources = "pubmed", enrich_compounds = FALSE,
    detail = "none", cache = FALSE,
    request_fun = function(url) stop("HTTP 503 service unavailable"),
    throttle = 0, progress = FALSE
  )
  expect_equal(busy$ProviderDiagnostics$rate_limit_count, 1L)
  expect_equal(busy$ProviderDiagnostics$timeout_count, 0L)
  expect_equal(busy$ProviderDiagnostics$no_hit_reason,
               "provider_error_or_incomplete_query")
  expect_match(busy$ProviderDiagnostics$warning_message, "503")
})

test_that("batch discovery preserves verified aliases and source identities", {
  root = tempfile("alias_batch_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  plants = data.frame(
    species = "Camellia sinensis",
    family = "Theaceae",
    verified_synonym = "Thea sinensis",
    stringsAsFactors = FALSE
  )
  lotus_index = data.frame(
    species = "Thea sinensis",
    genus = "Thea",
    family = "Theaceae",
    compound_name = "caffeine",
    lotus_id = "LTS_BATCH_ALIAS",
    cid = "2519",
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    inchikey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    molecular_formula = "C8H10N4O2",
    stringsAsFactors = FALSE
  )
  result = runPlantPhytochemistryBatch(
    plants = plants,
    sources = "lotus",
    taxon_fallback = "species",
    out_dir = root,
    cache_dir = file.path(root, "cache"),
    species_chunk_size = 1,
    compound_resolution_profile = "none",
    lotus_index = lotus_index,
    provider_indexes = list(lotus = lotus_index),
    require_all_providers = TRUE,
    resume = TRUE,
    overwrite = TRUE,
    progress = FALSE
  )
  expect_equal(nrow(result$PlantCompoundOccurrences), 1)
  expect_equal(result$PlantCompoundOccurrences$matched_taxon,
               "Thea sinensis")
  expect_true(any(result$PlantQueryAliases$alias == "Thea sinensis" &
                    result$PlantQueryAliases$verified == "Yes"))
  expect_equal(result$SourceCompoundIdentity$CID, 2519)
  expect_equal(result$BatchChunkManifest$status, "completed")

  aliases_changed = plants
  aliases_changed$verified_synonym = "Camellia thea"
  q1 = .plant_queries(plants, "species")
  q2 = .plant_queries(aliases_changed, "species")
  expect_false(identical(
    .plant_batch_query_signature(q1, .plant_query_aliases(plants, q1)),
    .plant_batch_query_signature(
      q2, .plant_query_aliases(aliases_changed, q2)
    )
  ))
})

test_that("multi-chunk batches retain one query identity per species", {
  root = tempfile("multichunk_query_ids_")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  plants = data.frame(
    species = paste("Simulata species", sprintf("%02d", seq_len(30))),
    verified_synonym = paste("Testophyta synonym", sprintf("%02d", seq_len(30))),
    stringsAsFactors = FALSE
  )
  literature = data.frame(
    species = plants$species,
    source_database = "PubMed",
    source_record_id = paste0("PMID", seq_len(30)),
    pmid = as.character(seq_len(30)),
    title = paste("Fixture record", seq_len(30)),
    evidence_tier = "direct_species_literature",
    confidence = "medium",
    stringsAsFactors = FALSE
  )

  result = runPlantPhytochemistryBatch(
    plants = plants,
    sources = "pubmed",
    out_dir = root,
    cache_dir = file.path(root, "cache"),
    species_chunk_size = 7,
    compound_resolution_profile = "none",
    provider_results = list(pubmed = list(
      LiteratureCandidates = literature
    )),
    resume = TRUE,
    overwrite = TRUE,
    progress = FALSE
  )

  accounting = result$ProviderQueryAccounting
  expect_equal(nrow(accounting), 30)
  expect_equal(length(unique(accounting$query_id)), 30)
  expect_setequal(accounting$species, plants$species)
  canonical_aliases = result$PlantQueryAliases[
    result$PlantQueryAliases$alias_type == "submitted_name", , drop = FALSE
  ]
  expect_equal(nrow(canonical_aliases), 30)
  expect_equal(length(unique(canonical_aliases$query_id)), 30)
  expect_setequal(canonical_aliases$species, plants$species)
})

test_that("identity resolution consumes non-LOTUS source structures", {
  occurrence = .plant_normalize_occurrences(data.frame(
    species = "Camellia sinensis",
    compound_name = "caffeine",
    source_database = "NPASS",
    source_record_id = "REC1",
    compound_id = "NPC0001",
    compound_id_type = "NPASS_NP_ID",
    evidence_tier = "direct_species_database",
    confidence = "high",
    stringsAsFactors = FALSE
  ))
  source_identity = data.frame(
    compound_name = "caffeine",
    compound_name_clean = "caffeine",
    source_database = "NPASS",
    source_record_id = "REC1",
    source_compound_id = "NPC0001",
    source_compound_id_type = "NPASS_NP_ID",
    CID = 2519L,
    InChIKey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N",
    SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    MolecularFormula = "C8H10N4O2",
    evidence_url = "https://bidd.group/NPASS/",
    evidence_text = "NPASS source record",
    identity_status = "source_structure_available",
    identity_note = "Fixture source structure.",
    stringsAsFactors = FALSE
  )
  input = list(
    PlantCompoundOccurrences = occurrence,
    SourceCompoundIdentity = source_identity
  )
  resolved = resolvePlantCompoundIdentities(
    input, source_only = TRUE,
    pubchem_fun = function(...) stop("PubChem must not be called"),
    progress = FALSE
  )
  expect_true(resolved$resolved)
  expect_equal(resolved$CID, "2519")
  expect_equal(resolved$InChIKey, "RYYVLZVUVIJVGH-UHFFFAOYSA-N")
  expect_match(resolved$resolution_source, "NPASS")
})

test_that("NCBI Datasets taxonomy suggestions require an exact verified taxon", {
  suggestions = list(sci_name_and_ids = list(
    list(sci_name = "Zoysia japonica", tax_id = "309978",
         matched_term = "Zoysia japonica", rank = "SPECIES"),
    list(sci_name = "Camellia japonica", tax_id = "4443",
         matched_term = "Camellia japonica", rank = "SPECIES")
  ))

  matched = .plant_verified_taxonomy_suggest_matches(
    suggestions, "Zoysia japonica"
  )
  absent = .plant_verified_taxonomy_suggest_matches(
    suggestions, "Geranium sessiliflorum"
  )

  expect_equal(nrow(matched), 1)
  expect_equal(matched$taxid, "309978")
  expect_equal(matched$scientific_name, "Zoysia japonica")
  expect_equal(nrow(absent), 0)
})

test_that("transient remote payloads are not accepted as valid caches", {
  cache_dir = tempfile("invalid_remote_cache_")
  on.exit(unlink(cache_dir, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(cache_dir, recursive = TRUE)
  url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed"
  cache_file = file.path(cache_dir, paste0(.pubchem_url_hash(url), ".json"))
  writeLines(
    '{"esearchresult":{"ERROR":"Search Backend failed: Status: 500"}}',
    cache_file
  )
  requests = 0L
  result = .plant_fetch_json(
    url, cache = TRUE, cache_dir = cache_dir, throttle = 0,
    request_fun = function(url) {
      requests <<- requests + 1L
      '{"esearchresult":{"count":"0","idlist":[]}}'
    }
  )

  expect_equal(requests, 1L)
  expect_false(.plant_cache_hit(result))
  expect_equal(.plant_pubmed_total_hits(result), 0L)
  expect_true(is.na(.plant_json_payload_error(result)))
  expect_true(.plant_service_busy_message("HTTP 500 returned by NCBI"))
  expect_match(
    .plant_text_payload_error(
      "<html><h1>Server Error</h1><div>Error: 500</div></html>"
    ),
    "500"
  )
})
