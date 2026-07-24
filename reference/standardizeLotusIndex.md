# Standardize a local LOTUS species-compound index

Converts a flat LOTUS export into a compact species-compound index that
uafR can query locally. This is the recommended LOTUS route for medium
or large plant panels because the live LOTUS simple API is unpaged and
can return very large payloads for common species.

The input may be a data frame or a path to a CSV, TSV, JSON, JSONL,
NDJSON, or RDS file produced from a LOTUS MongoDB, Wikidata, SDF
metadata, or other flat export. Column names are matched flexibly;
useful fields include \`species\`, \`genus\`, \`family\`, \`allTaxa\`,
\`traditional_name\`, \`iupac_name\`, \`lotus_id\`, \`inchikey\`,
\`smiles\`, \`molecular_formula\`, \`doi\`, \`pmid\`, \`plant_part\`,
\`tissue\`, and \`method\`.

## Usage

``` r
standardizeLotusIndex(x, source_file = NULL)
```

## Arguments

- x:

  Data frame or path to a flat LOTUS export.

- source_file:

  Optional source label recorded in the standardized index.

## Value

A normalized LOTUS index data frame with one row per taxon-compound
record where taxon evidence can be extracted.
