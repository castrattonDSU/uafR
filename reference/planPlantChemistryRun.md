# Plan a large plant chemistry run

Estimates run size, request burden, cache availability, pairwise output
sizes, and recommended conservative defaults before running a large
plant chemistry workflow. It does not query public services.

## Usage

``` r
planPlantChemistryRun(
  plants,
  compounds = NULL,
  sources = c("lotus", "pubmed", "pubtator"),
  cache_dir = NULL,
  lotus_index = NULL,
  expected_compounds_per_plant = 25,
  write_full_pairwise = FALSE,
  species_chunk_size = 25,
  compound_batch_size = 25,
  max_pubmed_records = 25
)
```

## Arguments

- plants:

  Character vector, data frame, or CSV path of plant names.

- compounds:

  Optional compound table/vector for categorate/Tanimoto estimates.

- sources:

  Provider sources expected for plant discovery.

- cache_dir:

  Optional cache directory to inspect.

- lotus_index:

  Optional local LOTUS index path.

- expected_compounds_per_plant:

  Estimated compounds per plant when no compound table is supplied.

- write_full_pairwise:

  Logical. If \`TRUE\`, estimate full pairwise file burden.

- species_chunk_size:

  Planned species discovery chunk size.

- compound_batch_size:

  Planned compound identity/enrichment batch size.

- max_pubmed_records:

  Planned maximum PubMed records per species.

## Value

List with \`Summary\`, \`PlantQueries\`, \`InputNameAudit\`,
\`ProviderPlan\`, \`CacheSummary\`, \`OutputEstimates\`,
\`ReadinessChecks\`, \`RunConfiguration\`, and \`Recommendations\`.
