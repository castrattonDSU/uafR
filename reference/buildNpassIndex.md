# Build a local NPASS species-metabolite index

Builds a manifest-backed, sharded lookup index from the official NPASS
3.0 natural-product general-information, structure, species-source, and
species taxonomy TSV files. The large species-source table is streamed
in chunks. No records or chemical identifiers are inferred: rows without
a source organism or named compound remain excluded from the lookup and
are counted in the build summary.

## Usage

``` r
buildNpassIndex(
  general_info,
  structures,
  species_pairs,
  species_info,
  out_dir,
  overwrite = FALSE,
  chunk_size = 1e+05,
  progress = interactive()
)
```

## Arguments

- general_info:

  Path to \`NPASS3.0_naturalproducts_generalinfo.txt\`.

- structures:

  Path to \`NPASS3.0_naturalproducts_structure.txt\`.

- species_pairs:

  Path to \`NPASS3.0_naturalproducts_species_pair.txt\`.

- species_info:

  Path to \`NPASS3.0_species_info.txt\`.

- out_dir:

  Output directory for sharded RDS files and \`manifest.json\`.

- overwrite:

  Logical. Replace an existing NPASS index directory.

- chunk_size:

  Number of species-pair lines parsed per streaming chunk.

- progress:

  Logical. Print chunk progress.

## Value

A list with \`BuildSummary\`, \`SourceManifest\`, \`ShardManifest\`, and
\`IndexPath\`.
