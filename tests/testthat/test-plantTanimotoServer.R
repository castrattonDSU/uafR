server_runner_path = function() {
  candidates = c(
    file.path(getwd(), "tools", "run_plant_tanimoto_server.R"),
    file.path(getwd(), "..", "tools", "run_plant_tanimoto_server.R"),
    testthat::test_path("..", "..", "tools",
                        "run_plant_tanimoto_server.R")
  )
  candidates = normalizePath(candidates, mustWork = FALSE)
  hit = candidates[file.exists(candidates)]
  if (length(hit) < 1L) return(NA_character_)
  hit[[1L]]
}

test_that("server runner completes all modes from cache without network", {
  runner = server_runner_path()
  skip_if(is.na(runner), "Server runner is not available in this test layout.")

  old_source = Sys.getenv("UAFR_SERVER_RUNNER_SOURCE_ONLY", unset = NA_character_)
  old_confirm = Sys.getenv("UAFR_CONFIRM_SERVER_TANIMOTO", unset = NA_character_)
  on.exit({
    if (is.na(old_source)) Sys.unsetenv("UAFR_SERVER_RUNNER_SOURCE_ONLY") else
      Sys.setenv(UAFR_SERVER_RUNNER_SOURCE_ONLY = old_source)
    if (is.na(old_confirm)) Sys.unsetenv("UAFR_CONFIRM_SERVER_TANIMOTO") else
      Sys.setenv(UAFR_CONFIRM_SERVER_TANIMOTO = old_confirm)
  }, add = TRUE)
  Sys.setenv(UAFR_SERVER_RUNNER_SOURCE_ONLY = "true",
             UAFR_CONFIRM_SERVER_TANIMOTO = "YES")
  runner_env = new.env(parent = globalenv())
  sys.source(runner, envir = runner_env)

  root = tempfile("server-runner-offline-")
  dir.create(root, recursive = TRUE)
  input_file = file.path(root, "membership.csv")
  species_universe_file = file.path(root, "plant_species_universe.csv")
  release_file = file.path(root, "uafR_release_manifest.json")
  out_dir = file.path(root, "output")
  cache_dir = file.path(root, "cache")
  input = data.frame(
    species = c("Plant A", "Plant B", "Plant C"),
    compound_id = c("cid_101", "cid_102", "cid_103"),
    compound_name = c("alpha", "beta", "gamma"),
    InChIKey = NA_character_,
    CID = c(101, 102, 103),
    comparison_scope = c("specialized_metabolites",
                         "specialized_metabolites", "primary_metabolites"),
    comparison_group = c("terpenoids", "terpenoids", "carbohydrates"),
    comparable_for_matrix = "Yes",
    stringsAsFactors = FALSE
  )
  utils::write.csv(input, input_file, row.names = FALSE, na = "")
  utils::write.csv(
    data.frame(species = c("Plant A", "Plant B", "Plant C", "Plant D")),
    species_universe_file, row.names = FALSE
  )
  jsonlite::write_json(
    list(package_version = as.character(utils::packageVersion("uafR"))),
    release_file, auto_unbox = TRUE
  )
  fingerprint = function(character) paste(rep(character, 154), collapse = "")

  url = paste0(
    .pubchem_base_url(),
    "/pug/compound/cid/101,102,103/property/",
    paste(c("Fingerprint2D", "CanonicalSMILES", "IsomericSMILES",
            "InChIKey", "MolecularFormula", "Title"), collapse = ","),
    "/JSON"
  )
  properties = list(PropertyTable = list(Properties = list(
    list(CID = 101, Title = "alpha", MolecularFormula = "CH4",
         InChIKey = "AAAAAAAAAAAAAA-BBBBBBBBBB-C",
         CanonicalSMILES = "C", IsomericSMILES = "C",
         Fingerprint2D = fingerprint("B")),
    list(CID = 102, Title = "beta", MolecularFormula = "C2H6",
         InChIKey = "CCCCCCCCCCCCCC-DDDDDDDDDD-E",
         CanonicalSMILES = "CC", IsomericSMILES = "CC",
         Fingerprint2D = fingerprint("C")),
    list(CID = 103, Title = "gamma", MolecularFormula = "C3H8",
         InChIKey = "FFFFFFFFFFFFFF-GGGGGGGGGG-H",
         CanonicalSMILES = "CCC", IsomericSMILES = "CCC",
         Fingerprint2D = fingerprint("/"))
  )))
  tanimoto_cache = file.path(cache_dir, "tanimoto")
  dir.create(tanimoto_cache, recursive = TRUE)
  stem = file.path(tanimoto_cache, .pubchem_url_hash(url))
  jsonlite::write_json(properties, paste0(stem, ".json"), auto_unbox = TRUE)
  writeLines(url, paste0(stem, ".url"), useBytes = TRUE)

  base_args = c("--input", input_file, "--out-dir", out_dir,
                "--cache-dir", cache_dir,
                "--species-universe", species_universe_file,
                "--release-manifest", release_file,
                "--throttle", "1.1")
  expect_identical(runner_env$main(c(base_args, "--mode", "preflight")), 0L)
  expect_identical(runner_env$main(c(base_args, "--mode", "smoke",
                                     "--smoke-structures", "3")), 0L)
  expect_identical(runner_env$main(c(base_args, "--mode", "summary")), 0L)
  expect_identical(runner_env$main(c(base_args, "--mode", "full",
                                     "--full-pairs", "true")), 0L)

  expect_true(file.exists(file.path(out_dir,
                                    "FINGERPRINT_RESOLUTION_COMPLETED.txt")))
  expect_true(file.exists(file.path(out_dir, "SUMMARY_COMPLETED.txt")))
  expect_true(file.exists(file.path(out_dir, "FULL_PAIRWISE_COMPLETED.txt")))
  expect_false(file.exists(file.path(out_dir, "FAILED.json")))
  expect_false(file.exists(file.path(out_dir, "PAUSED_SERVICE_BUSY.json")))
  expect_false(any(grepl("[.]partial$", list.files(out_dir, recursive = TRUE))))
  manifest = utils::read.csv(
    file.path(out_dir, "server_tanimoto_output_manifest.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_setequal(unique(manifest$mode), c("smoke", "summary", "full"))
  validation = utils::read.csv(
    file.path(out_dir, "full", "server_tanimoto_validation.csv"),
    stringsAsFactors = FALSE
  )
  expect_true(all(validation$status == "pass"))
  complete_summary = utils::read.csv(
    file.path(out_dir, "full", "plant_pair_tanimoto_summary.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(nrow(complete_summary), choose(4, 2))
  expect_equal(sum(complete_summary$support_status == "insufficient_support"),
               3L)
  expect_true(file.exists(file.path(
    out_dir, "full", "comparable_scope_tanimoto_summary.csv"
  )))
  scope_summary = utils::read.csv(
    file.path(out_dir, "full", "comparable_scope_tanimoto_summary.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expect_equal(scope_summary$comparison_scope, "specialized_metabolites")
  preflight_checks = utils::read.csv(
    file.path(out_dir, "server_tanimoto_preflight_checks.csv"),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  release_check = preflight_checks[
    preflight_checks$check == "installed_package_matches_release_manifest",
    , drop = FALSE
  ]
  expect_equal(nrow(release_check), 1L)
  expect_identical(release_check$status, "pass")

  bad_release = file.path(root, "bad_release_manifest.json")
  bad_out = file.path(root, "bad_release_output")
  jsonlite::write_json(list(package_version = "0.0.0"), bad_release,
                       auto_unbox = TRUE)
  expect_error(
    runner_env$main(c("--input", input_file, "--out-dir", bad_out,
                      "--cache-dir", cache_dir,
                      "--release-manifest", bad_release,
                      "--mode", "preflight", "--throttle", "1.1")),
    "preflight failed"
  )
  expect_false(file.exists(file.path(bad_out, "server_progress.csv")))

  busy_out = file.path(root, "service_busy_output")
  busy_cache = file.path(root, "service_busy_cache")
  busy_args = c("--input", input_file, "--out-dir", busy_out,
                "--cache-dir", busy_cache, "--throttle", "1.1",
                "--mode", "summary", "--service-busy-limit", "1")
  busy_status = runner_env$main(
    busy_args,
    request_fun = function(url) stop("HTTP status 503")
  )
  expect_identical(busy_status, 75L)
  expect_true(file.exists(file.path(busy_out, "PAUSED_SERVICE_BUSY.json")))
  expect_false(file.exists(file.path(busy_out, "FAILED.json")))
  busy_state = jsonlite::read_json(
    file.path(busy_out, "server_run_status.json"), simplifyVector = TRUE
  )
  expect_identical(busy_state$state, "paused_service_busy")
  archive_root = file.path(busy_out, ".incomplete", "archive")
  expect_true(dir.exists(archive_root))
  expect_true(length(list.dirs(archive_root, recursive = FALSE)) > 0L)
})
