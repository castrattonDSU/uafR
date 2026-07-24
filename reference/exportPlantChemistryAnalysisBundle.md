# Export a plant chemistry analysis bundle from batch outputs

\`exportPlantChemistryAnalysisBundle()\` writes a researcher-facing
bundle that combines resumable categorate batches with optional
plant-chemistry and Tanimoto tables. It is designed for large
species-first projects where enrichment batches, species-compound
memberships, and pairwise chemistry summaries need to be handed to a
downstream analysis workspace.

## Usage

``` r
exportPlantChemistryAnalysisBundle(
  categorate_batches,
  path,
  plant_membership = NULL,
  species_pair_tanimoto = NULL,
  resolved_compounds = NULL,
  pubchem_fingerprints = NULL,
  plant_compound_pair_tanimoto = NULL,
  file_references = NULL,
  plant_list = NULL,
  metadata = NULL,
  project_id = NULL,
  tables = NULL,
  format = c("auto", "xlsx", "csv"),
  include_empty = FALSE,
  min_property_ratio = 0.9,
  overwrite = FALSE,
  max_cell_chars = 30000,
  finalize = TRUE,
  include_feature_exports = TRUE,
  include_comparable_tanimoto = FALSE,
  validate_export = TRUE
)
```

## Arguments

- categorate_batches:

  A batch directory path, object returned by
  \`readCategorateBatchDirectory()\`, list of categorate result objects,
  or character vector of \`.rds\` files.

- path:

  Output directory for \`format = "csv"\` or \`.xlsx\` path for \`format
  = "xlsx"\`.

- plant_membership:

  Optional data frame or CSV path containing species-compound membership
  rows.

- species_pair_tanimoto:

  Optional data frame or CSV path containing species-pair Tanimoto
  summaries.

- resolved_compounds:

  Optional data frame or CSV path containing unique resolved compounds.

- pubchem_fingerprints:

  Optional data frame or CSV path containing PubChem fingerprint
  records.

- plant_compound_pair_tanimoto:

  Optional data frame or CSV/CSV.GZ path containing plant-compound pair
  Tanimoto rows. This can be large, so it is read only when comparable
  Tanimoto summaries are requested.

- file_references:

  Optional character vector or data frame of large external files to
  record in the bundle manifest, such as compressed compound-pair
  Tanimoto tables.

- plant_list:

  Optional character vector, data frame, or CSV path containing the full
  plant list used by the project. When supplied, the finalized bundle
  includes a missing-chemistry coverage table.

- metadata:

  Optional plant metadata data frame or CSV path. If it contains
  \`species\`, \`accepted_species_name\`, \`genus\`, or \`family\`,
  those fields are used only for transparent taxonomy/status joins and
  are never fabricated.

- project_id:

  Optional project label written to bundle documentation and provenance
  text.

- tables:

  Categorate result table names to combine and export. If \`NULL\`, a
  curated analysis-ready set is used.

- format:

  Export format: \`"csv"\`, \`"xlsx"\`, or \`"auto"\`.

- include_empty:

  Logical. If \`TRUE\`, export empty combined tables.

- min_property_ratio:

  Minimum acceptable ratio of \`PubChemProperties\` rows to resolved
  PubChem CIDs.

- overwrite:

  Logical. If \`TRUE\`, replace an existing output.

- max_cell_chars:

  Maximum characters retained in a single character cell.

- finalize:

  Logical. If \`TRUE\` and \`format = "csv"\`, add enriched
  analysis-ready tables, bundle documentation, and validation summaries.

- include_feature_exports:

  Logical. If \`TRUE\`, finalized CSV bundles include model-ready
  species feature matrices.

- include_comparable_tanimoto:

  Logical. If \`TRUE\`, finalized CSV bundles include scope- and
  group-filtered Tanimoto summaries. Requires
  \`plant_compound_pair_tanimoto\`.

- validate_export:

  Logical. If \`TRUE\`, run bundle CSV validation after finalization.

## Value

A manifest data frame describing exported tables.

## Examples

``` r
if (FALSE) { # \dontrun{
manifest = exportPlantChemistryAnalysisBundle(
  categorate_batches = "categorate_batches",
  path = "plant_chemistry_analysis_bundle",
  plant_membership = "dsi_species_compound_membership.csv",
  species_pair_tanimoto = "dsi_species_pair_tanimoto_summary.csv",
  overwrite = TRUE
)
manifest
} # }
```
