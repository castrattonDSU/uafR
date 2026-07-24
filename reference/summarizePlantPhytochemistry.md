# Summarize plant phytochemistry

Summarize plant phytochemistry

## Usage

``` r
summarizePlantPhytochemistry(
  plant_compounds,
  categorate_result = NULL,
  compound_resolution = NULL,
  plant_queries = NULL,
  provider_diagnostics = NULL,
  comparability = NULL
)
```

## Arguments

- plant_compounds:

  Plant-compound occurrence table or a \`uaf_plant_phytochemistry\`
  result.

- categorate_result:

  Optional categorate-like enrichment result.

- compound_resolution:

  Optional compound resolution table.

- plant_queries:

  Optional plant query table.

- provider_diagnostics:

  Optional provider diagnostics table.

- comparability:

  Optional precomputed \`ChemistryComparability\` table.

## Value

Species-level summary table.
