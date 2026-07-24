# Build a compact local LOTUS index

Builds the compact LOTUS index used by \`queryLotusIndex()\` and the
\`lotus_index\` argument in plant phytochemistry workflows. Input can be
one or more flat LOTUS exports, a directory containing supported flat
exports, or an in-memory data frame. Supported file formats are CSV,
TSV, JSON, JSONL, NDJSON, and RDS.

Raw LOTUS MongoDB downloads are commonly distributed as zipped BSON
dumps. Use \`tools/flatten_lotus_mongo_dump.py\` to stream the official
LOTUS MongoDB ZIP into an auditable flat CSV, a compact CSV, and
optionally a manifest-backed lookup directory. Then pass the
flat/compact file to this function or pass the lookup directory directly
to \`queryLotusIndex()\` or \`resolvePlantPhytochemistry(lotus_index =
...)\`.

## Usage

``` r
buildLotusIndex(
  input,
  out_file = NULL,
  format = c("auto", "csv", "rds"),
  overwrite = FALSE,
  manifest_file = NULL,
  recursive = TRUE,
  strict = FALSE,
  progress = interactive()
)
```

## Arguments

- input:

  Data frame, file path, directory path, or character vector of
  file/directory paths containing flat LOTUS export rows.

- out_file:

  Optional output path for the compact LOTUS index.

- format:

  Output format. \`"auto"\` infers from \`out_file\`; otherwise use
  \`"csv"\` or \`"rds"\`.

- overwrite:

  Logical. If \`FALSE\`, existing output files are not replaced.

- manifest_file:

  Optional JSON manifest path. If omitted and \`out_file\` is supplied,
  a sidecar \`\*\_manifest.json\` file is written.

- recursive:

  Logical. If \`TRUE\`, directory inputs are searched recursively for
  supported flat export files.

- strict:

  Logical. If \`TRUE\`, stop on the first unreadable input source. If
  \`FALSE\`, unreadable sources are recorded in the build manifest and
  other sources are still processed.

- progress:

  Logical. If \`TRUE\`, print source-processing progress.

## Value

A list with \`LotusIndex\`, \`BuildSummary\`, and \`BuildManifest\`.
