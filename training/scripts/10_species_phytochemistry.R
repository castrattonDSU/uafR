# Species-first plant phytochemistry smoke test.
#
# This script is offline by default. It demonstrates the curated intake path,
# validation, species summaries, matrices, candidate scoring, and export. Live
# PubMed discovery and compound enrichment are optional because they depend on
# internet access and public-service availability.

run_species_phytochemistry_smoke <- function(output_dir = NULL,
                                             cache_dir = NULL) {
  if (!requireNamespace("uafR", quietly = TRUE)) {
    if (file.exists("DESCRIPTION") && requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(".", quiet = TRUE)
    } else {
      stop("uafR is not installed. Install uafR or run this script from the ",
           "uafR source tree with devtools available.", call. = FALSE)
    }
  } else {
    library(uafR)
  }

  plants <- c("Salix nigra", "Camellia sinensis", "Zea mays")
  curated <- data.frame(
    species = plants,
    compound_name = c("salicin", "caffeine", "DIMBOA"),
    source_database = c("training_curated", "training_curated",
                        "training_curated"),
    citation_or_url = c("training example", "training example",
                        "training example"),
    evidence_tier = c("manual_curated", "manual_curated", "manual_curated"),
    plant_part = c("bark", "leaf", "seedling"),
    method = c("LC-MS", "GC-MS", "LC-MS"),
    evidence_note = c(
      "Simulated training row for the curated intake workflow.",
      "Simulated training row for the curated intake workflow.",
      "Simulated training row for the curated intake workflow."
    ),
    stringsAsFactors = FALSE
  )

  phyto <- resolvePlantPhytochemistry(
    plants = plants,
    sources = character(),
    curated_data = curated,
    enrich_compounds = FALSE,
    detail = "none"
  )

  validation <- validatePlantPhytochemistryResult(phyto)
  analysis_ready <- filterPlantPhytochemistryEvidence(phyto)
  matrix <- plantPhytochemistryMatrix(phyto, mode = "binary")
  scores <- scorePlantChemistryCandidates(phyto)

  message("Offline smoke test: unresolved-compound warnings are expected ",
          "because live compound enrichment is disabled.")
  print(validation$Summary)
  print(phyto$SpeciesChemistrySummary)
  print(analysis_ready$PlantCompoundOccurrences)
  print(matrix)
  print(scores)

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    exportPlantPhytochemistryWorkbook(
      phyto,
      path = file.path(output_dir, "plant_phytochemistry_export"),
      format = "csv",
      overwrite = TRUE
    )
    exportPlantPhytochemistryWorkbook(
      phyto,
      path = file.path(output_dir, "plant_phytochemistry_analysis_ready_export"),
      format = "csv",
      preset = "analysis_ready",
      overwrite = TRUE
    )
  }

  if (identical(tolower(Sys.getenv("UAFR_LIVE_PLANT_DISCOVERY")), "true")) {
    if (is.null(cache_dir)) {
      cache_dir <- file.path(tempdir(), "uafR_plant_cache")
    }
    live <- resolvePlantPhytochemistry(
      plants = plants,
      sources = c("pubmed"),
      enrich_compounds = FALSE,
      detail = "none",
      cache = TRUE,
      cache_dir = cache_dir,
      max_pubmed_records = 5
    )
    print(live$ProviderDiagnostics)
    print(utils::head(live$LiteratureCandidates))
  }

  invisible(phyto)
}

if (identical(environment(), globalenv())) {
  run_species_phytochemistry_smoke()
}
