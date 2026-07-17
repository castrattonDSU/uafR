fixture_pubchem_request = function(url) {
  if (grepl("/pug/compound/name/aspirin/cids/JSON$", url)) {
    return(list(IdentifierList = list(CID = list(2244))))
  }

  if (grepl("/property/", url)) {
    return(list(PropertyTable = list(Properties = list(list(
      CID = 2244,
      Title = "Aspirin",
      MolecularFormula = "C9H8O4",
      MolecularWeight = "180.159",
      IUPACName = "2-acetyloxybenzoic acid",
      InChIKey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
      CanonicalSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O",
      IsomericSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O",
      ExactMass = "180.04225873",
      TPSA = "63.6",
      Complexity = "212"
    )))))
  }

  if (grepl("/synonyms/JSON$", url)) {
    return(list(InformationList = list(Information = list(list(
      CID = 2244,
      Synonym = list("aspirin", "acetylsalicylic acid")
    )))))
  }

  if (grepl("heading=Mass%20Spectrometry|heading=GC-MS|heading=MS-MS", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 1,
        SourceName = "Fixture Spectra",
        URL = "https://example.test/spectra"
      )),
      Section = list(list(
        TOCHeading = "Mass Spectrometry",
        Section = list(list(
          TOCHeading = "GC-MS",
          Information = list(list(
            Name = "Top m/z peaks",
            ReferenceNumber = 1,
            Value = list(StringWithMarkup = list(
              list(String = "120 100"),
              list(String = "92 80")
            ))
          ))
        ))
      ))
    )))
  }

  if (grepl("heading=Safety%20and%20Hazards", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 2,
        SourceName = "Fixture Safety",
        URL = "https://example.test/safety"
      )),
      Section = list(list(
        TOCHeading = "Safety and Hazards",
        Section = list(list(
          TOCHeading = "GHS Classification",
          Information = list(
            list(
              Name = "Signal",
              ReferenceNumber = 2,
              Value = list(StringWithMarkup = list(list(String = "Warning")))
            ),
            list(
              Name = "Hazard Statement",
              ReferenceNumber = 2,
              Value = list(StringWithMarkup = list(list(String = "Causes serious eye irritation")))
            )
          )
        ))
      ))
    )))
  }

  if (grepl("heading=Experimental%20Properties", url)) {
    return(list(Record = list(
      Reference = list(list(
        ReferenceNumber = 3,
        SourceName = "Fixture Experimental",
        URL = "https://example.test/experimental"
      )),
      Section = list(list(
        TOCHeading = "Experimental Properties",
        Information = list(list(
          Name = "Boiling Point",
          ReferenceNumber = 3,
          Value = list(StringWithMarkup = list(list(String = "140 C")))
        ))
      ))
    )))
  }

  if (grepl("/assaysummary/JSON$", url)) {
    return(list(Table = list(Row = list(list(AID = 123, ActivityOutcome = "Active")))))
  }

  if (grepl("/assay/aid/123/description/JSON$", url)) {
    return(list(PC_AssayContainer = list(list(assay = list(descr = list(
      aid = list(id = 123, version = 1),
      aid_source = list(db = list(name = "Fixture BioAssay",
                                  source_id = list(str = "FIX123"))),
      name = "qHTS Assay for Inhibitors of PTGS1",
      description = list("Assay Overview: PTGS1 inhibition assay"),
      protocol = list("Measure prostaglandin synthesis inhibition."),
      comment = list("Keywords: inhibitor, cyclooxygenase"),
      xref = list(
        list(xref = list(pmid = 123456)),
        list(xref = list(gene = 5743)),
        list(xref = list(taxonomy = 9606)),
        list(xref = list(dburl = "https://example.test/bioassay"))
      ),
      results = list(
        list(tid = 1, name = "Potency"),
        list(tid = 2, name = "Efficacy")
      ),
      target = list(list(
        name = "prostaglandin G/H synthase 1 [Homo sapiens]",
        descr = "Human PTGS1 (COX-1)",
        mol_id = list(protein_accession = "NP_000953.2"),
        organism = list(org = list(
          taxname = "Homo sapiens",
          common = "human",
          db = list(list(db = "taxon", tag = list(id = 9606)))
        ))
      )),
      activity_outcome_method = 2,
      project_category = 1
    ))))))
  }

  list(Record = list(Section = list()))
}

test_that("PUG-View information without a usable value has a stable row schema", {
  information = list(list(
    Name = "Unavailable annotation",
    ReferenceNumber = 1,
    Value = list(StringWithMarkup = list(list(String = "")))
  ))
  references = list(`1` = list(
    Source = "Fixture source",
    URL = "https://example.test/unavailable"
  ))

  rows = uafR:::.pubchem_parse_information(
    information = information,
    cid = 2244,
    query = "aspirin",
    heading = "Chemical and Physical Properties",
    heading_path = "Chemical and Physical Properties > Fixture",
    references = references,
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/2244"
  )

  expect_length(rows, 1)
  expect_s3_class(rows[[1]], "data.frame")
  expect_equal(nrow(rows[[1]]), 1)
  expect_true(all(c("Value", "MarkupText", "MarkupURL", "MarkupExtra") %in%
                    names(rows[[1]])))
  expect_true(all(is.na(rows[[1]][1, c(
    "Value", "MarkupText", "MarkupURL", "MarkupExtra"
  )])))
})

test_that("PubChem cache keys are collision-resistant and request-bound", {
  inchikeys = c(
    "DFYRUELUNQRZTB-UHFFFAOYSA-N",
    "OILXMJHPFNGGTO-ZAUYPBDWSA-N",
    "VTZLLYJSIAQXFJ-UHFFFAOYSA-N"
  )
  urls = paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/",
    inchikeys, "/cids/JSON"
  )
  hashes = vapply(urls, .pubchem_url_hash, character(1))
  expect_length(unique(hashes), length(urls))
  expect_true(all(grepl("^v2_[0-9a-f]{32}$", hashes)))

  calls = new.env(parent = emptyenv())
  calls$n = 0L
  request_fun = function(url) {
    calls$n = calls$n + 1L
    cid = match(url, urls) + 100L
    jsonlite::toJSON(list(IdentifierList = list(CID = list(cid))),
                     auto_unbox = TRUE)
  }
  cache_dir = tempfile("pubchem-cache-v2-")
  fetch = .pubchem_fetcher(cache = TRUE, cache_dir = cache_dir,
                           throttle = 0, request_fun = request_fun)
  first = lapply(urls, fetch)
  second = lapply(urls, fetch)

  expect_equal(calls$n, length(urls))
  expect_equal(vapply(first, .tanimoto_json_cid, integer(1)), 101:103)
  expect_equal(vapply(second, .tanimoto_json_cid, integer(1)), 101:103)
  expect_equal(length(list.files(cache_dir, pattern = "[.]json$")), 3L)
  expect_equal(length(list.files(cache_dir, pattern = "[.]url$")), 3L)
  bindings = vapply(
    paste0(file.path(cache_dir, hashes), ".url"),
    function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
    character(1)
  )
  expect_identical(unname(bindings), urls)
})

test_that("pubchemProfile returns structured profile tables from PubChem responses", {
  profile = pubchemProfile("aspirin",
                           profile = "full",
                           cache = FALSE,
                           throttle = 0,
                           request_fun = fixture_pubchem_request)

  expect_s3_class(profile, "uaf_pubchem_profile")
  expect_equal(profile$identity$CID, 2244)
  expect_equal(profile$properties$MolecularFormula, "C9H8O4")
  expect_true("TPSA" %in% colnames(profile$properties))
  expect_true("acetylsalicylic acid" %in% profile$synonyms$Synonym)
  expect_true(any(profile$spectra$Value == "120 100"))
  expect_true(any(profile$safety$Value == "Warning"))
  expect_true(any(profile$experimental$Name == "Boiling Point"))
  expect_true("CleanValue" %in% colnames(profile$annotations))
  expect_equal(profile$experimental$ValueNumeric[profile$experimental$Name == "Boiling Point"][[1]],
               140)
  expect_true(any(profile$bioactivity$Value == "Active"))
  expect_true(all(c("AID", "ActivityOutcome", "ActivityClass",
                    "ActivityValue", "AssayName", "TargetAccession") %in%
                    colnames(profile$bioactivity)))
  expect_true(any(profile$bioactivity$ActivityClass == "active"))
  expect_equal(nrow(profile$bioassay_details), 1)
  expect_true(any(profile$bioassay_details$TargetName ==
                    "prostaglandin G/H synthase 1 [Homo sapiens]" &
                    profile$bioassay_details$TargetGeneID == "5743" &
                    profile$bioassay_details$TargetTaxonomyID == "9606"))
  expect_true(nrow(profile$provenance) > 0)
})

test_that("pubchemProfile resolves InChIKey queries through the InChIKey endpoint", {
  inchikey = "GWCVYMLSGTZAHQ-HNNXBMFYSA-N"
  request_fun = function(url) {
    if (grepl("/pug/compound/name/", url)) {
      stop("InChIKey query should not use the name endpoint")
    }
    if (grepl("/pug/compound/inchikey/GWCVYMLSGTZAHQ-HNNXBMFYSA-N/cids/JSON$",
              url)) {
      return(list(IdentifierList = list(CID = list(162905822))))
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 162905822,
        Title = "Fixture InChIKey compound",
        MolecularFormula = "C22H23NO6",
        InChIKey = inchikey,
        CanonicalSMILES = "COc1ccccc1",
        IsomericSMILES = "COc1ccccc1"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 162905822,
        Synonym = list("fixture compound")
      )))))
    }
    list()
  }

  expect_no_warning(
    profile <- pubchemProfile(
      inchikey,
      profile = "minimal",
      include_annotations = FALSE,
      request_fun = request_fun
    )
  )

  expect_equal(profile$identity$CID, 162905822)
  expect_equal(profile$identity$MatchStatus, "resolved_inchikey")
  expect_match(profile$identity$SourceURL, "/compound/inchikey/")
  expect_equal(profile$properties$InChIKey, inchikey)
})

test_that("pubchemProfile accepts explicit PubChem CID query tokens", {
  requested = character()
  request_fun = function(url) {
    requested <<- c(requested, url)
    if (grepl("/pug/compound/name/", url) ||
        grepl("/pug/compound/inchikey/", url)) {
      stop("CID query should not use name or InChIKey endpoints")
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 2244,
        Title = "Aspirin",
        MolecularFormula = "C9H8O4",
        InChIKey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        CanonicalSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O",
        IsomericSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 2244,
        Synonym = list("aspirin")
      )))))
    }
    list()
  }

  expect_no_warning(
    profile <- pubchemProfile(
      "cid:2244",
      profile = "minimal",
      include_annotations = FALSE,
      request_fun = request_fun
    )
  )

  expect_equal(profile$identity$CID, 2244)
  expect_equal(profile$identity$MatchStatus, "resolved_cid")
  expect_match(profile$identity$SourceURL, "/compound/cid/2244/cids/JSON")
  expect_equal(profile$properties$Title, "Aspirin")
})

test_that("pubchemProfile falls back to single-CID property requests", {
  request_fun = function(url) {
    if (grepl("/property/", url) &&
        grepl("/cid/2244,2519/property/", url)) {
      return(list())
    }
    if (grepl("/property/", url) &&
        grepl("/cid/2244/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 2244,
        Title = "Aspirin",
        MolecularFormula = "C9H8O4",
        InChIKey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
      )))))
    }
    if (grepl("/property/", url) &&
        grepl("/cid/2519/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 2519,
        Title = "Caffeine",
        MolecularFormula = "C8H10N4O2",
        InChIKey = "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      cid = if (grepl("/cid/2244/", url)) 2244 else 2519
      return(list(InformationList = list(Information = list(list(
        CID = cid,
        Synonym = list(paste0("synonym-", cid))
      )))))
    }
    list()
  }

  profile = pubchemProfile(
    c("cid:2244", "cid:2519"),
    profile = "minimal",
    include_annotations = FALSE,
    request_fun = request_fun
  )

  expect_equal(sort(profile$properties$CID), c(2244, 2519))
  expect_true(all(c("Aspirin", "Caffeine") %in% profile$properties$Title))
})

test_that("pubchemProfile validates empty compound input", {
  expect_error(pubchemProfile(c("", NA), request_fun = fixture_pubchem_request),
               "at least one non-empty")
})

test_that("pubchemProfile can skip PUG-View annotations for property-only workflows", {
  property_only_request = function(url) {
    if (grepl("/pug_view/", url)) {
      stop("PUG-View should not be requested when include_annotations = FALSE")
    }
    fixture_pubchem_request(url)
  }

  profile = pubchemProfile("aspirin",
                           profile = "ms",
                           cache = FALSE,
                           throttle = 0,
                           include_annotations = FALSE,
                           request_fun = property_only_request)

  expect_equal(profile$identity$CID, 2244)
  expect_equal(profile$properties$MolecularFormula, "C9H8O4")
  expect_equal(nrow(profile$annotations), 0)
  expect_equal(nrow(profile$source_annotations), 0)
})

test_that("pubchemProfile can fetch and locally filter one full PUG-View record", {
  requested = character()
  record_request = function(url) {
    requested <<- c(requested, url)
    if (grepl("/pug/compound/name/aspirin/cids/JSON$", url)) {
      return(list(IdentifierList = list(CID = list(2244))))
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 2244, Title = "Aspirin", MolecularFormula = "C9H8O4",
        InChIKey = "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
        CanonicalSMILES = "CC(=O)OC1=CC=CC=C1C(=O)O"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 2244, Synonym = list("aspirin")
      )))))
    }
    if (grepl("/pug_view/data/compound/2244/JSON$", url)) {
      return(list(Record = list(
        Reference = list(
          list(ReferenceNumber = 1, SourceName = "Fixture Safety",
               URL = "https://example.test/safety"),
          list(ReferenceNumber = 2,
               SourceName = "LOTUS - the natural products occurrence database",
               URL = "https://example.test/lotus")
        ),
        Section = list(
          list(
            TOCHeading = "Safety and Hazards",
            Information = list(list(
              Name = "Signal", ReferenceNumber = 1,
              Value = list(StringWithMarkup = list(list(String = "Warning")))
            ))
          ),
          list(
            TOCHeading = "Natural Products",
            Information = list(list(
              Name = "Occurrence", ReferenceNumber = 2,
              Value = list(StringWithMarkup = list(list(
                String = "Reported natural-product record"
              )))
            ))
          )
        )
      )))
    }
    stop("Unexpected request: ", url)
  }

  profile = pubchemProfile(
    "aspirin", profile = "minimal", sections = "Safety and Hazards",
    sources = "LOTUS - the natural products occurrence database",
    annotation_request_mode = "record", cache = FALSE, throttle = 0,
    request_fun = record_request
  )

  pug_view = grepl("/pug_view/data/compound/2244/JSON", requested,
                   fixed = TRUE)
  expect_equal(sum(pug_view), 1L)
  expect_false(any(grepl("?heading=", requested, fixed = TRUE)))
  expect_false(any(grepl("?source=", requested, fixed = TRUE)))
  expect_true(any(profile$annotations$HeadingPath == "Safety and Hazards"))
  expect_true(any(grepl("LOTUS", profile$source_annotations$Source,
                        fixed = TRUE)))
  expect_identical(profile$annotation_request_mode, "record")
})

test_that("record parsing materializes only requested headings and sources", {
  json = list(Record = list(
    Reference = list(
      list(ReferenceNumber = 1, SourceName = "Fixture Safety",
           URL = "https://example.test/safety"),
      list(ReferenceNumber = 2,
           SourceName = "LOTUS - the natural products occurrence database",
           URL = "https://example.test/lotus"),
      list(ReferenceNumber = 3, SourceName = "Fixture Literature",
           URL = "https://example.test/literature")
    ),
    Section = list(
      list(
        TOCHeading = "Safety and Hazards",
        Information = list(list(
          Name = "Signal", ReferenceNumber = 1,
          Value = list(StringWithMarkup = list(list(String = "Warning")))
        ))
      ),
      list(
        TOCHeading = "Natural Products",
        Information = list(list(
          Name = "Occurrence", ReferenceNumber = 2,
          Value = list(StringWithMarkup = list(list(
            String = "Reported natural-product record"
          )))
        ))
      ),
      list(
        TOCHeading = "Literature",
        Information = list(list(
          Name = "Reference", ReferenceNumber = 3,
          Value = list(StringWithMarkup = list(list(
            String = "Irrelevant full-record literature row"
          )))
        ))
      )
    )
  ))

  full = uafR:::.pubchem_parse_pugview(
    json, cid = 2244, query = "aspirin",
    heading = "Full PubChem record",
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/2244"
  )
  selected = uafR:::.pubchem_parse_pugview(
    json, cid = 2244, query = "aspirin",
    heading = "Full PubChem record",
    pubchem_url = "https://pubchem.ncbi.nlm.nih.gov/compound/2244",
    include_headings = "Safety and Hazards",
    include_sources = "LOTUS - the natural products occurrence database"
  )
  expected = uafR:::.pubchem_bind_tables(
    uafR:::.pubchem_select_annotation_headings(
      full, "Safety and Hazards"
    ),
    uafR:::.pubchem_select_annotation_sources(
      full, "LOTUS - the natural products occurrence database"
    )
  )

  expect_equal(selected, expected)
  expect_setequal(selected$HeadingPath,
                  c("Safety and Hazards", "Natural Products"))
  expect_false(any(grepl("Irrelevant", selected$CleanValue, fixed = TRUE)))
})

test_that("classification parsing retains large source-backed hierarchy sets", {
  hierarchy = function(index) {
    list(
      SourceName = "LOTUS - the natural products occurrence database",
      HID = 115,
      Information = list(HID = 115),
      Node = list(
        list(
          NodeID = paste0("leaf_", index),
          ParentID = paste0("parent_", index),
          Information = list(
            HNID = paste0("hnid_", index),
            Match = TRUE,
            Name = list(StringWithMarkup = list(list(
              String = paste("Species", index)
            )))
          )
        ),
        list(
          NodeID = paste0("parent_", index),
          Information = list(
            HNID = paste0("parent_hnid_", index),
            Name = list(StringWithMarkup = list(list(String = "Plantae")))
          )
        )
      )
    )
  }
  json = list(Hierarchies = list(
    Hierarchy = lapply(seq_len(1000), hierarchy)
  ))
  link = data.frame(
    Query = "fixture compound", CID = 1L,
    Source = "LOTUS - the natural products occurrence database",
    TreeID = "115", TreeName = "Biological Classification",
    TreeType = "biological", SourceURL = "https://example.test/source",
    PubChemURL = "https://example.test/compound",
    stringsAsFactors = FALSE
  )

  parsed = uafR:::.pubchem_parse_classification_response(
    json, link, "https://example.test/classification"
  )

  expect_equal(nrow(parsed), 1000L)
  expect_equal(length(unique(parsed$ClassName)), 1000L)
  expect_true(all(parsed$ParentClass == "Plantae"))
  expect_true(all(parsed$ClassDepth == 2L))
})

test_that("pubchemProfile strict requests preserve retryable failures", {
  expect_error(
    pubchemProfile(
      "aspirin", profile = "minimal", cache = FALSE, throttle = 0,
      service_busy_limit = 2, max_attempts = 3,
      fail_on_retry_exhausted = TRUE,
      request_fun = function(url) stop("HTTP status 503")
    ),
    class = "uaf_pubchem_service_busy"
  )
})

test_that("pubchemProfile can omit synonyms from identity-only requests", {
  no_synonym_request = function(url) {
    if (grepl("/synonyms/JSON$", url)) stop("Synonyms should be skipped")
    fixture_pubchem_request(url)
  }
  profile = pubchemProfile(
    "aspirin", profile = "minimal", include_annotations = FALSE,
    include_synonyms = FALSE, cache = FALSE, throttle = 0,
    request_fun = no_synonym_request
  )
  expect_equal(profile$identity$CID, 2244L)
  expect_equal(profile$properties$MolecularFormula, "C9H8O4")
  expect_equal(nrow(profile$synonyms), 0L)
})

test_that("pubchemProfile retries conservative name aliases for PubChem identity", {
  requested = character()
  alias_request = function(url) {
    requested <<- c(requested, utils::URLdecode(url))
    if (grepl("/pug/compound/name/beta-pinene/cids/JSON$", utils::URLdecode(url))) {
      return(list(IdentifierList = list(CID = list(14896))))
    }
    if (grepl("/property/", url)) {
      return(list(PropertyTable = list(Properties = list(list(
        CID = 14896,
        Title = "beta-Pinene",
        MolecularFormula = "C10H16",
        MolecularWeight = "136.238",
        IUPACName = "6,6-dimethyl-2-methylidenebicyclo[3.1.1]heptane",
        InChIKey = "WTARULDDTDQWMU-UHFFFAOYSA-N",
        SMILES = "CC1(C2CCC(=C)C1C2)C"
      )))))
    }
    if (grepl("/synonyms/JSON$", url)) {
      return(list(InformationList = list(Information = list(list(
        CID = 14896,
        Synonym = list("beta-pinene")
      )))))
    }
    list()
  }
  compound = paste0(intToUtf8(0x03b2), "-pinene")

  profile = pubchemProfile(compound,
                           profile = "minimal",
                           cache = FALSE,
                           throttle = 0,
                           include_annotations = FALSE,
                           request_fun = alias_request)

  expect_equal(profile$identity$CID, 14896)
  expect_equal(profile$identity$MatchStatus, "resolved_alias")
  expect_equal(profile$identity$QueriedName, "beta-pinene")
  expect_true(any(grepl("/name/beta-pinene/cids/JSON", requested,
                        fixed = TRUE)))
  expect_equal(profile$properties$SMILES, "CC1(C2CCC(=C)C1C2)C")
})

test_that("PubChem fetcher opens and resets its service-busy circuit breaker", {
  old_verbose = Sys.getenv("UAFR_PUBCHEM_VERBOSE", unset = NA_character_)
  on.exit({
    if (is.na(old_verbose)) {
      Sys.unsetenv("UAFR_PUBCHEM_VERBOSE")
    } else {
      Sys.setenv(UAFR_PUBCHEM_VERBOSE = old_verbose)
    }
  }, add = TRUE)
  Sys.setenv(UAFR_PUBCHEM_VERBOSE = "false")

  events = list()
  always_busy = .pubchem_fetcher(
    cache = FALSE,
    cache_dir = tempfile("pubchem-busy-"),
    throttle = 0,
    request_fun = function(url) stop("HTTP status 503"),
    service_busy_limit = 2,
    event_fun = function(event) events[[length(events) + 1L]] <<- event,
    max_attempts = 3
  )
  expect_error(always_busy("https://example.test/busy"),
               class = "uaf_pubchem_service_busy")
  event_names = vapply(events, `[[`, character(1), "event")
  expect_equal(sum(event_names == "request_attempt"), 2L)
  expect_true("service_busy_limit_reached" %in% event_names)

  attempt = 0L
  busy_then_success = .pubchem_fetcher(
    cache = FALSE,
    cache_dir = tempfile("pubchem-reset-"),
    throttle = 0,
    request_fun = function(url) {
      attempt <<- attempt + 1L
      if (attempt %% 2L == 1L) stop("HTTP status 503")
      list(value = attempt)
    },
    service_busy_limit = 2,
    max_attempts = 2
  )
  expect_equal(busy_then_success("https://example.test/one")$value, 2L)
  expect_equal(busy_then_success("https://example.test/two")$value, 4L)
})

test_that("PubChem fetcher treats 404 as a no-hit rather than service failure", {
  events = list()
  fetch = .pubchem_fetcher(
    cache = FALSE,
    cache_dir = tempfile("pubchem-no-hit-"),
    throttle = 0,
    request_fun = function(url) stop("HTTP status 404"),
    service_busy_limit = 1,
    event_fun = function(event) events[[length(events) + 1L]] <<- event,
    max_attempts = 1
  )
  expect_null(fetch("https://example.test/not-found"))
  event_names = vapply(events, `[[`, character(1), "event")
  expect_true("request_no_hit" %in% event_names)
  expect_false("service_busy_limit_reached" %in% event_names)
})
