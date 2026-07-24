# Export model-ready plant chemistry feature tables

\`exportPlantChemistryFeatureSet()\` builds analysis-neutral feature
tables from enriched plant-compound membership rows. The tables are
suitable for joining to phylogenetic, ecological, water-quality,
remediation, or other project metadata.

## Usage

``` r
exportPlantChemistryFeatureSet(
  membership,
  path = NULL,
  overwrite = FALSE,
  modes = "count",
  include_context = TRUE,
  include_evidence = TRUE,
  species_universe = NULL,
  plant_metadata = NULL
)
```

## Arguments

- membership:

  Enriched plant-compound membership table.

- path:

  Optional output directory. If supplied, CSV files and a manifest are
  written.

- overwrite:

  Logical. If \`TRUE\`, replace an existing output directory.

- modes:

  Matrix modes to export. Supported values are \`"count"\`,
  \`"binary"\`, \`"fraction"\`, and \`"confidence"\`. The default
  preserves the original count-matrix output.

- include_context:

  Logical. If \`TRUE\`, include plant-part, tissue, and method matrices
  when those fields are present.

- include_evidence:

  Logical. If \`TRUE\`, include evidence-grade matrices.

- species_universe:

  Optional character vector or data frame defining the complete species
  set and row order for every feature matrix. Species present in
  \`membership\` but absent from this input are appended. The default
  uses all species in \`membership\`.

- plant_metadata:

  Optional plant metadata table used to populate taxonomy fields for
  species with no membership rows.

## Value

Named list of feature tables and manifest.
