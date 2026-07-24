# Query a local NPASS species-metabolite index

Query a local NPASS species-metabolite index

## Usage

``` r
queryNpassIndex(
  plants,
  npass_index,
  taxon_fallback = c("species", "genus"),
  max_records = Inf
)
```

## Arguments

- plants:

  Character vector or data frame of plant names.

- npass_index:

  Manifest-backed directory created by \`buildNpassIndex()\`, or a
  standardized NPASS index data frame.

- taxon_fallback:

  Taxonomic ranks to query. Species is always queried; genus and family
  rows remain explicit fallback evidence.

- max_records:

  Maximum records retained per input plant.

## Value

A normalized \`PlantCompoundOccurrences\` table. Source structure
evidence is attached as the \`SourceCompoundIdentity\` attribute.
