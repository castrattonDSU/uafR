fixture_kegg_request = function(url) {
  compound_record = paste(c(
    "ENTRY       C01405                      Compound",
    "NAME        Aspirin;",
    "            Acetylsalicylic acid",
    "FORMULA     C9H8O4",
    "EXACT_MASS  180.0423",
    "MOL_WEIGHT  180.159",
    "PATHWAY     map00590  Arachidonic acid metabolism",
    "REACTION    R01335",
    "ENZYME      3.1.1.55",
    "DBLINKS     PubChem: 2244",
    "            ChEBI: 15365"
  ), collapse = "\n")

  drug_record = paste(c(
    "ENTRY       D00109                      Drug",
    "NAME        Aspirin (JP18/USP);",
    "            Acetylsalicylic acid",
    "FORMULA     C9H8O4",
    "MOL_WEIGHT  180.159",
    "BRITE       Anatomical Therapeutic Chemical (ATC) classification",
    "DBLINKS     PubChem: 2244"
  ), collapse = "\n")

  if (grepl("/find/compound/aspirin$", url)) {
    return("cpd:C01405\tAspirin; Acetylsalicylic acid")
  }
  if (grepl("/find/drug/aspirin$", url)) {
    return("dr:D00109\tAspirin (JP18/USP)")
  }
  if (grepl("/get/C01405\\+D00109$", url) ||
      grepl("/get/D00109\\+C01405$", url)) {
    return(paste(compound_record, "///", drug_record, "///", sep = "\n"))
  }
  if (grepl("/get/C01405$", url)) return(paste(compound_record, "///", sep = "\n"))
  if (grepl("/get/D00109$", url)) return(paste(drug_record, "///", sep = "\n"))
  if (grepl("/get/path:map00590$", url)) {
    return(paste(c(
      "ENTRY       map00590                    Pathway",
      "NAME        Arachidonic acid metabolism",
      "///"
    ), collapse = "\n"))
  }
  if (grepl("/get/rn:R01335$", url)) {
    return(paste(c(
      "ENTRY       R01335                      Reaction",
      "NAME        acetylsalicylate deacetylase reaction",
      "DEFINITION  Aspirin + H2O <=> Salicylate + Acetate",
      "EQUATION    C01405 + C00001 <=> C00805 + C00033",
      "///"
    ), collapse = "\n"))
  }
  if (grepl("/get/ec:3.1.1.55$", url)) {
    return(paste(c(
      "ENTRY       EC 3.1.1.55                 Enzyme",
      "NAME        acetylsalicylate deacetylase",
      "///"
    ), collapse = "\n"))
  }

  if (grepl("/link/pathway/cpd:C01405$", url)) {
    return("cpd:C01405\tpath:map00590")
  }
  if (grepl("/link/reaction/cpd:C01405$", url)) {
    return("cpd:C01405\trn:R01335")
  }
  if (grepl("/link/enzyme/cpd:C01405$", url)) {
    return("cpd:C01405\tec:3.1.1.55")
  }
  if (grepl("/link/pubmed/cpd:C01405$", url)) {
    return("cpd:C01405\tpubmed:123456")
  }
  ""
}

test_that("keggProfile parses KEGG records, links, identifiers, and groups", {
  profile = keggProfile("aspirin",
                        kegg_ids = "KEGG: C01405",
                        cache = FALSE,
                        throttle = 0,
                        request_fun = fixture_kegg_request)

  expect_s3_class(profile, "uaf_kegg_profile")
  expect_true("C01405" %in% profile$matches$KEGG_ID)
  expect_true(any(profile$records$Field == "FORMULA" &
                    profile$records$CleanValue == "C9H8O4"))
  expect_true(any(profile$identifiers$IdentifierType == "PubChem" &
                    profile$identifiers$Identifier == "2244"))
  expect_true(any(profile$pathways$PathwayID == "map00590"))
  expect_true(any(profile$pathways$PathwayName == "Arachidonic acid metabolism"))
  expect_true(any(profile$pathways$PathwayGroup == "Lipid metabolism",
                  na.rm = TRUE))
  expect_true(any(profile$reactions$ReactionID == "R01335"))
  expect_true(any(profile$reactions$ReactionName == "acetylsalicylate deacetylase reaction",
                  na.rm = TRUE))
  expect_true(any(profile$reactions$Equation == "C01405 + C00001 <=> C00805 + C00033",
                  na.rm = TRUE))
  expect_true(any(profile$reactions$ReactionDefinition == "Aspirin + H2O <=> Salicylate + Acetate",
                  na.rm = TRUE))
  expect_true(any(profile$enzymes$EnzymeName == "acetylsalicylate deacetylase",
                  na.rm = TRUE))
  expect_true(any(profile$enzymes$EnzymeClass == "Hydrolases",
                  na.rm = TRUE))
  expect_true(any(profile$link_metadata$TargetDatabase == "reaction" &
                    profile$link_metadata$Name == "acetylsalicylate deacetylase reaction"))
  expect_true(any(profile$classifications$Classification == "Lipid metabolism"))
  expect_true(nrow(profile$provenance) > 0)
})

test_that("keggProfile validates empty input", {
  expect_error(keggProfile(cache = FALSE,
                           throttle = 0,
                           request_fun = fixture_kegg_request),
               "Provide at least one")
})

test_that("keggProfile uses KEGG-safe keyword paths for punctuated names", {
  requested = character()
  fixture = function(url) {
    requested <<- c(requested, url)
    if (grepl("/find/compound/1\\+2\\+4-Trimethylbenzene$", url)) {
      return("cpd:C14533\t1,2,4-Trimethylbenzene; Pseudocumene")
    }
    ""
  }

  profile = keggProfile(
    "1,2,4-Trimethylbenzene", cache = FALSE, throttle = 0,
    link_targets = character(), request_fun = fixture
  )

  expect_true(any(profile$matches$KEGG_ID == "C14533"))
  expect_true(any(grepl("/find/compound/1\\+2\\+4-Trimethylbenzene$",
                        requested)))
  expect_false(any(grepl("%2C", requested, fixed = TRUE)))
  expect_equal(.kegg_encode_find_query("alpha beta"), "alpha+beta")
})

test_that("KEGG broad substring hits remain rejected search candidates", {
  fixture = function(url) {
    if (grepl("/find/compound/1-Hexanol$", url)) {
      return(paste(
        "cpd:C02498\t2-Ethylhexan-1-ol; 2-Ethyl-1-hexanol",
        "cpd:C00854\t1-Hexanol; Hexan-1-ol",
        sep = "\n"
      ))
    }
    ""
  }

  profile = keggProfile(
    "1-Hexanol", cache = FALSE, throttle = 0,
    link_targets = character(), request_fun = fixture
  )

  expect_equal(profile$matches$KEGG_ID, "C00854")
  expect_equal(profile$matches$MatchStatus, "exact_name_match")
  rejected = profile$search_candidates[
    profile$search_candidates$KEGG_ID == "C02498", , drop = FALSE
  ]
  expect_equal(rejected$Accepted, "No")
  expect_equal(rejected$RejectionReason, "broad_name_match_not_exact")
  expect_true(is.na(rejected$MatchRank))
})

test_that("KEGG fetcher caches no-record 404 and propagates service failures", {
  cache_dir = tempfile("kegg_http_cache_")
  on.exit(unlink(cache_dir, recursive = TRUE, force = TRUE), add = TRUE)
  calls = 0L
  no_record = function(url) {
    calls <<- calls + 1L
    stop("HTTP 404 returned for ", url, call. = FALSE)
  }
  fetch = .kegg_fetcher(
    cache = TRUE, cache_dir = cache_dir, throttle = 0,
    request_fun = no_record
  )

  expect_identical(fetch("https://rest.kegg.jp/get/br:br08011"), "")
  expect_identical(fetch("https://rest.kegg.jp/get/br:br08011"), "")
  expect_equal(calls, 1L)

  busy = .kegg_fetcher(
    cache = FALSE, cache_dir = cache_dir, throttle = 0,
    request_fun = function(url) {
      stop("HTTP 503 returned for ", url, call. = FALSE)
    }
  )
  expect_error(busy("https://rest.kegg.jp/info/kegg"), "HTTP 503")
  expect_equal(.kegg_http_status("HTTP status was '429 Too Many Requests'"),
               429L)
})

test_that("keggProfile limits broad name-search expansion", {
  requested = new.env(parent = emptyenv())
  requested$urls = character()

  broad_fixture = function(url) {
    requested$urls = c(requested$urls, url)

    if (grepl("/find/compound/glucose$", url)) {
      return(paste(c(
        "cpd:C00031\tD-Glucose; Glucose",
        "cpd:C00267\talpha-D-Glucose; Glucose",
        "cpd:C00221\tbeta-D-Glucose; Glucose",
        "cpd:C01172\tD-Glucose 6-phosphate",
        "cpd:C00668\talpha-D-Glucose 6-phosphate"
      ), collapse = "\n"))
    }
    if (grepl("/find/drug/glucose$", url)) return("")

    if (grepl("/get/", url)) {
      ids = strsplit(sub(".*/get/", "", url), "\\+")[[1]]
      records = vapply(ids, function(id) {
        paste(c(
          paste0("ENTRY       ", id, "                      Compound"),
          paste0("NAME        Fixture ", id),
          "FORMULA     C6H12O6",
          "MOL_WEIGHT  180.156"
        ), collapse = "\n")
      }, character(1))
      return(paste(paste(records, collapse = "\n///\n"), "///", sep = "\n"))
    }

    if (grepl("/link/pathway/", url)) {
      id = sub(".*/cpd:", "", url)
      return(paste0("cpd:", id, "\tpath:map00010"))
    }

    ""
  }

  profile = keggProfile("glucose",
                        cache = FALSE,
                        throttle = 0,
                        max_matches_per_query = 2,
                        link_targets = "pathway",
                        request_fun = broad_fixture)

  expect_lte(sum(profile$matches$Query == "glucose"), 2)
  expect_equal(sort(profile$matches$MatchRank), c(1L, 2L))
  expect_true(all(profile$links$TargetDatabase == "pathway"))
  expect_false(any(grepl("/link/enzyme/", requested$urls)))
  expect_lte(sum(grepl("/link/pathway/", requested$urls)), 2)
})
