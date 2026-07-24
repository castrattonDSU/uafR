# Inspect plant-provider readiness

Reports whether each requested species-first provider has the local
index or live-service configuration required for a run. This is a
configuration and resource check by default; it does not make live
requests.

## Usage

``` r
plantProviderAvailability(
  sources = c("lotus", "knapsack", "npass", "pubchem", "pubmed", "pubtator"),
  provider_indexes = NULL,
  lotus_index = Sys.getenv("UAFR_LOTUS_INDEX", ""),
  probe_live = FALSE,
  request_fun = NULL,
  request_timeout = 15
)
```

## Arguments

- sources:

  Provider names.

- provider_indexes:

  Named list of local provider indexes. Supported names currently
  include \`lotus\` and \`npass\`.

- lotus_index:

  Backward-compatible LOTUS index argument.

- probe_live:

  Logical. If \`TRUE\`, perform a lightweight request to live provider
  landing/API endpoints.

- request_fun:

  Optional request function used when \`probe_live = TRUE\`.

- request_timeout:

  Maximum live probe time in seconds.

## Value

A provider resource manifest data frame.
