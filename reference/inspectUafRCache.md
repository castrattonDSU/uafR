# Inspect a uafR cache directory

Lists cache files and summarizes size, extension, and likely provider
from path names. This function never queries live services.

## Usage

``` r
inspectUafRCache(cache_dir, recursive = TRUE)
```

## Arguments

- cache_dir:

  Cache directory.

- recursive:

  Logical. If \`TRUE\`, inspect nested files.

## Value

Cache inventory data frame.
