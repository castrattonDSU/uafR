# Finalize a plant chemistry analysis bundle

\`finalizePlantChemistryAnalysisBundle()\` adds the analysis-ready layer
to a CSV bundle written by \`exportPlantChemistryAnalysisBundle()\`. It
does not query public services. Instead, it uses the exported
plant-compound membership, Tanimoto summaries, PubChem
fingerprints/properties, and categorate-derived tables already present
in the bundle.

The finalizer writes enriched plant-compound membership, species
chemistry summaries, missing-species coverage, source coverage
summaries, validation overviews, a data dictionary, and README/methods
text. It also refreshes the export manifest so downstream projects can
audit row and column counts.

## Usage

``` r
finalizePlantChemistryAnalysisBundle(
  path,
  plant_list = NULL,
  metadata = NULL,
  plant_compound_pair_tanimoto = NULL,
  project_id = NULL,
  overwrite = TRUE,
  include_feature_exports = TRUE,
  include_comparable_tanimoto = FALSE,
  validate_export = TRUE,
  max_cell_chars = 30000
)
```

## Arguments

- path:

  CSV bundle directory.

- plant_list:

  Optional character vector, data frame, or CSV path containing the
  complete project plant list. A \`species\` column is preferred when a
  data frame is supplied.

- metadata:

  Optional plant metadata data frame or CSV path with a \`species\`
  column. \`accepted_species_name\`, \`genus\`, and \`family\` are
  joined when available.

- plant_compound_pair_tanimoto:

  Optional plant-compound pair Tanimoto data frame or CSV/CSV.GZ path
  used to build comparable scope/group summaries when
  \`include_comparable_tanimoto = TRUE\`.

- project_id:

  Optional project label written to documentation.

- overwrite:

  Logical. If \`FALSE\`, existing finalization outputs are not
  overwritten.

- include_feature_exports:

  Logical. If \`TRUE\`, write model-ready species feature matrices and
  species quality metadata.

- include_comparable_tanimoto:

  Logical. If \`TRUE\`, write Tanimoto summaries filtered to matched
  comparable chemistry scopes and groups.

- validate_export:

  Logical. If \`TRUE\`, run bundle validation and include validation
  status in \`12b_ValidationOverview.csv\`.

- max_cell_chars:

  Maximum characters retained in a single exported cell.

## Value

Updated export manifest.
