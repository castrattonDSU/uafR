representative_pilot_tool_path = function() {
  candidates = c(
    file.path(getwd(), "tools", "run_representative_lotus_pilot.R"),
    file.path(getwd(), "..", "tools", "run_representative_lotus_pilot.R"),
    testthat::test_path("..", "..", "tools",
                        "run_representative_lotus_pilot.R")
  )
  candidates = normalizePath(candidates, mustWork = FALSE)
  hit = candidates[file.exists(candidates)]
  if (length(hit) < 1L) return(NA_character_)
  hit[[1L]]
}

test_that("representative LOTUS pilot selection is deterministic and stratified", {
  tool = representative_pilot_tool_path()
  skip_if(is.na(tool), "Representative pilot tool is unavailable.")
  old = Sys.getenv("UAFR_REPRESENTATIVE_PILOT_SOURCE_ONLY",
                   unset = NA_character_)
  on.exit({
    if (is.na(old)) {
      Sys.unsetenv("UAFR_REPRESENTATIVE_PILOT_SOURCE_ONLY")
    } else {
      Sys.setenv(UAFR_REPRESENTATIVE_PILOT_SOURCE_ONLY = old)
    }
  }, add = TRUE)
  Sys.setenv(UAFR_REPRESENTATIVE_PILOT_SOURCE_ONLY = "true")
  env = new.env(parent = globalenv())
  sys.source(tool, envir = env)

  plants = data.frame(
    species = sprintf("Plant species %02d", 1:36),
    lotus_exact_species_key_available = c(rep(TRUE, 12), rep(FALSE, 24)),
    lotus_genus_key_available = c(rep(TRUE, 24), rep(FALSE, 12)),
    stringsAsFactors = FALSE
  )
  panel_a = env$select_panel(
    plants, "species", "lotus_exact_species_key_available",
    "lotus_genus_key_available", 10, 10, 5
  )
  panel_b = env$select_panel(
    plants[36:1, ], "species", "lotus_exact_species_key_available",
    "lotus_genus_key_available", 10, 10, 5
  )

  expect_equal(nrow(panel_a), 25L)
  expect_equal(anyDuplicated(panel_a$species), 0L)
  expect_identical(panel_a$species, panel_b$species)
  counts = table(factor(
    panel_a$pilot_stratum,
    levels = c("exact_species_key", "genus_key_no_exact_species_key",
               "no_species_or_genus_key")
  ))
  expect_identical(as.integer(counts), c(10L, 10L, 5L))
})
