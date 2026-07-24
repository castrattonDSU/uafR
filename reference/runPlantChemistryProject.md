# Run a reusable plant chemistry project workflow

\`runPlantChemistryProject()\` coordinates plant list/metadata intake,
curated or discovered plant-compound rows, optional cached categorate
batch exports, bundle finalization, validation, evidence grading, and
feature-set outputs. Live discovery is opt-in through \`run_discovery =
TRUE\`.

## Usage

``` r
runPlantChemistryProject(
  plant_list,
  output_dir,
  metadata = NULL,
  plant_compounds = NULL,
  plant_result = NULL,
  categorate_batches = NULL,
  species_pair_tanimoto = NULL,
  resolved_compounds = NULL,
  pubchem_fingerprints = NULL,
  plant_compound_pair_tanimoto = NULL,
  file_references = NULL,
  project_id = NULL,
  run_discovery = FALSE,
  overwrite = FALSE,
  ...
)
```

## Arguments

- plant_list:

  Character vector, data frame, or CSV path of plant names.

- output_dir:

  Project output directory.

- metadata:

  Optional plant metadata table or path.

- plant_compounds:

  Optional curated plant-compound table or path.

- plant_result:

  Optional result object or RDS path from plant workflows.

- categorate_batches:

  Optional categorate batch directory or files.

- species_pair_tanimoto:

  Optional plant-pair Tanimoto summary table/path.

- resolved_compounds:

  Optional resolved compounds table/path.

- pubchem_fingerprints:

  Optional PubChem fingerprints table/path.

- plant_compound_pair_tanimoto:

  Optional plant-compound pair table/path for comparable Tanimoto
  summaries.

- file_references:

  Optional file references recorded in the bundle.

- project_id:

  Optional project label.

- run_discovery:

  Logical. If \`TRUE\`, call \`resolvePlantPhytochemistry()\`.

- overwrite:

  Logical. If \`TRUE\`, replace existing project output.

- ...:

  Additional arguments passed to \`resolvePlantPhytochemistry()\` when
  \`run_discovery = TRUE\`.

## Value

Project manifest list.
