# Validate a plant chemistry analysis bundle

\`validatePlantChemistryAnalysisBundle()\` checks a CSV bundle for
parser consistency, manifest row/column agreement, portable artifact
paths, file sizes and checksums, accidental row-index columns, and
required columns in the main analysis tables. Python \`csv.reader\` and
pandas checks are optional so package tests do not depend on a Python
installation, but the function records whether those checks were run.

## Usage

``` r
validatePlantChemistryAnalysisBundle(
  path,
  use_python = FALSE,
  use_pandas = FALSE
)
```

## Arguments

- path:

  CSV bundle directory.

- use_python:

  Logical. If \`TRUE\`, also validate every CSV with Python's standard
  \`csv.reader\` when \`python3\` is available.

- use_pandas:

  Logical. If \`TRUE\`, also validate every CSV with
  \`pandas.read_csv()\` when Python and pandas are available.

## Value

A list with \`Summary\`, \`CSVValidation\`, \`RequiredColumns\`,
\`ArtifactValidation\`, \`ManifestReferences\`, and \`Manifest\` tables.
