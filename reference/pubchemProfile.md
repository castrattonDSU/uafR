# Build a detailed PubChem profile for query chemicals

`pubchemProfile()` retrieves structured PubChem information for query
chemicals and returns tidy tables that can be reused by downstream uafR
workflows. It separates PubChem enrichment from mass spectrometry
processing so results can be cached, inspected, and tested
independently.

## Usage

``` r
pubchemProfile(
  compounds,
  profile = c("minimal", "ms", "safety", "bioactivity", "full"),
  sections = NULL,
  sources = NULL,
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.2,
  assay_detail_limit = 50,
  request_fun = NULL
)
```

## Arguments

- compounds:

  A character vector of chemical names to resolve on PubChem.

- profile:

  One of `"minimal"`, `"ms"`, `"safety"`, `"bioactivity"`, or `"full"`.
  Profiles control which PubChem properties and PUG-View headings are
  requested.

- sections:

  Optional PUG-View headings to request in addition to the selected
  profile.

- sources:

  Optional PUG-View source names to request in addition to the selected
  profile. Use this for source-specific records such as LOTUS, FEMA,
  FDA/SPL, or MeSH.

- cache:

  Logical. If `TRUE`, raw PubChem JSON responses are cached under
  `cache_dir`.

- cache_dir:

  Directory for cached PubChem responses. Defaults to a user-cache
  location from
  [`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html) when
  available, otherwise a temporary directory.

- throttle:

  Seconds to wait between uncached PubChem requests.

- assay_detail_limit:

  Maximum number of unique PubChem BioAssay AIDs per query for which
  assay-description metadata is fetched when `profile` is
  `"bioactivity"` or `"full"`. Active, numeric, and target-bearing
  assays are prioritized. Set to `0` to skip assay-description requests.

- request_fun:

  Optional function used to retrieve a URL. This is intended for tests
  and advanced users. It should return either JSON text or a parsed
  list.

## Value

A list with tidy data frames: `identity`, `properties`, `synonyms`,
`annotations`, `source_annotations`, `taxonomy`, `classifications`,
`spectra`, `safety`, `experimental`, `bioactivity`, `bioassay_details`,
and `provenance`. Annotation-derived tables include raw `Value`,
normalized `CleanValue`, parsed `ValueNumeric`, and normalized
`UnitClean` columns when those fields can be extracted. The
`bioactivity` table stores one PubChem BioAssay summary row per assay
with activity outcome, assay name/type, target identifiers, and numeric
activity values when PubChem reports them. The `bioassay_details` table
stores bounded assay-description metadata for selected AIDs. The
returned object has class `"uaf_pubchem_profile"`.

## Examples

``` r
if (FALSE) { # \dontrun{
profile = pubchemProfile(c("methyl salicylate", "octanal"), profile = "ms")
profile$properties
profile$spectra
} # }
```
