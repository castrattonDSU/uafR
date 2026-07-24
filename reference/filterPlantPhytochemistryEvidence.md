# Filter plant phytochemistry occurrence evidence

Keeps plant-compound occurrence rows that meet explicit evidence,
confidence, and biological-context criteria. When \`x\` is a plant
phytochemistry result, the returned result is rebuilt with updated
occurrence, compound-resolution, trait-evidence, summary, matrix, and
validation tables. Literature candidate tables are preserved for
audit/provenance.

## Usage

``` r
filterPlantPhytochemistryEvidence(
  x,
  occurrence_status = c("direct_reported", "curated_reported"),
  analysis_ready = TRUE,
  min_confidence = "medium",
  plant_part_group = NULL,
  tissue_group = NULL,
  method_group = NULL,
  source_database = NULL,
  evidence_tier = NULL
)
```

## Arguments

- x:

  Plant phytochemistry result or normalized occurrence table.

- occurrence_status:

  Character vector of occurrence status values to retain. Defaults to
  direct and curated reported records. Use \`"all"\` or \`NULL\` to skip
  this filter.

- analysis_ready:

  Optional logical. If \`TRUE\`, retain only rows marked analysis-ready;
  if \`FALSE\`, retain only non-analysis-ready rows; if \`NULL\`, do not
  filter on the flag.

- min_confidence:

  Minimum confidence to retain. Use \`NULL\` to skip.

- plant_part_group:

  Optional plant-part groups to retain.

- tissue_group:

  Optional tissue groups to retain.

- method_group:

  Optional method groups to retain.

- source_database:

  Optional source databases to retain.

- evidence_tier:

  Optional evidence tiers to retain.

## Value

A filtered plant phytochemistry result when \`x\` is a result object;
otherwise a filtered \`PlantCompoundOccurrences\` data frame.
