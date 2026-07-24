# Query a local LOTUS index for plant phytochemistry records

Searches a standardized local LOTUS index for direct species records
and, when requested, genus or family fallback records. The returned
table uses the same \`PlantCompoundOccurrences\` schema as
\`resolvePlantPhytochemistry()\`.

## Usage

``` r
queryLotusIndex(
  plants,
  lotus_index,
  taxon_fallback = c("species", "genus"),
  max_records = Inf
)
```

## Arguments

- plants:

  Character vector or data frame of plant names.

- lotus_index:

  Data frame or path accepted by \`standardizeLotusIndex()\`, or a
  manifest-backed lookup directory produced by
  \`tools/flatten_lotus_mongo_dump.py –lookup-dir\`.

- taxon_fallback:

  Fallback ranks. Species records are always queried; optional
  \`"genus"\` and \`"family"\` records are labeled as fallback evidence.

- max_records:

  Maximum records retained per input plant.

## Value

A normalized \`PlantCompoundOccurrences\` data frame.
