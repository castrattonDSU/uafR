# Build a detailed PubChem profile for query chemicals

\`pubchemProfile()\` retrieves structured PubChem information for query
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
  include_annotations = TRUE,
  include_synonyms = TRUE,
  annotation_request_mode = c("filtered", "record"),
  service_busy_limit = Inf,
  request_event_fun = NULL,
  max_attempts = NULL,
  fail_on_retry_exhausted = FALSE,
  query_overrides = NULL,
  request_fun = NULL
)
```

## Arguments

- compounds:

  A character vector of chemical names to resolve on PubChem.

- profile:

  One of \`"minimal"\`, \`"ms"\`, \`"safety"\`, \`"bioactivity"\`, or
  \`"full"\`. Profiles control which PubChem properties and PUG-View
  headings are requested.

- sections:

  Optional PUG-View headings to request in addition to the selected
  profile.

- sources:

  Optional PUG-View source names to request in addition to the selected
  profile. Use this for source-specific records such as LOTUS, FEMA,
  FDA/SPL, or MeSH.

- cache:

  Logical. If \`TRUE\`, raw PubChem JSON responses are cached under
  \`cache_dir\`.

- cache_dir:

  Directory for cached PubChem responses. Defaults to a user-cache
  location from \`tools::R_user_dir()\` when available, otherwise a
  temporary directory.

- throttle:

  Seconds to wait between uncached PubChem requests.

- assay_detail_limit:

  Maximum number of unique PubChem BioAssay AIDs per query for which
  assay-description metadata is fetched when \`profile\` is
  \`"bioactivity"\` or \`"full"\`. Active, numeric, and target-bearing
  assays are prioritized. Set to \`0\` to skip assay-description
  requests.

- include_annotations:

  Logical. If \`TRUE\`, fetch PubChem PUG-View annotation sections
  selected by \`profile\`, \`sections\`, and \`sources\`. Set to
  \`FALSE\` when only identity and property fields are needed.

- include_synonyms:

  Logical. If \`TRUE\`, fetch PubChem synonyms for resolved CIDs.
  Identity-only production runs can set this to \`FALSE\` to avoid an
  unnecessary request while retaining identity and property fields.

- annotation_request_mode:

  One of \`"filtered"\` or \`"record"\`. \`"filtered"\` makes one
  PUG-View request per requested heading and source. \`"record"\`
  retrieves one complete PUG-View record per CID and filters the parsed
  record locally. The latter is intended for large resumable runs and
  can substantially reduce request counts, at the cost of larger
  responses.

- service_busy_limit:

  Number of consecutive HTTP 429/503 responses allowed before a
  \`uaf_pubchem_service_busy\` condition stops the request. Use \`Inf\`
  to preserve the general interactive default.

- request_event_fun:

  Optional callback receiving one list per PubChem request/cache event.
  Intended for operational monitoring.

- max_attempts:

  Optional maximum attempts per uncached request.

- fail_on_retry_exhausted:

  Logical. If \`TRUE\`, exhausted retryable or transport failures raise
  a \`uaf_pubchem_request_failed\` condition instead of being treated as
  an ordinary no-hit. Production batch workflows should enable this so
  temporary service failures remain retryable.

- query_overrides:

  Optional named character vector or two-column data frame mapping each
  displayed compound name to a source-backed PubChem query, such as
  \`"cid:2519"\` or a full InChIKey. The displayed \`Query\` remains the
  compound name. Overrides do not fabricate identity and should come
  from a reviewed source record.

- request_fun:

  Optional function used to retrieve a URL. This is intended for tests
  and advanced users. It should return either JSON text or a parsed
  list.

## Value

A list with tidy data frames: \`identity\`, \`properties\`,
\`synonyms\`, \`annotations\`, \`source_annotations\`, \`taxonomy\`,
\`classifications\`, \`spectra\`, \`safety\`, \`experimental\`,
\`bioactivity\`, \`bioassay_details\`, and \`provenance\`.
Annotation-derived tables include raw \`Value\`, normalized
\`CleanValue\`, parsed \`ValueNumeric\`, and normalized \`UnitClean\`
columns when those fields can be extracted. The \`identity\` table
includes the exact \`QueriedName\` sent to PubChem and marks
conservative deterministic alias matches, such as Greek-letter
transliterations and trailing-punctuation cleanup, with \`MatchStatus =
"resolved_alias"\`. The \`bioactivity\` table stores one PubChem
BioAssay summary row per assay with activity outcome, assay name/type,
target identifiers, and numeric activity values when PubChem reports
them. The \`bioassay_details\` table stores bounded assay-description
metadata for selected AIDs. The returned object has class
\`"uaf_pubchem_profile"\`.

## Examples

``` r
if (FALSE) { # \dontrun{
profile = pubchemProfile(c("methyl salicylate", "octanal"), profile = "ms")
profile$properties
profile$spectra
} # }
```
