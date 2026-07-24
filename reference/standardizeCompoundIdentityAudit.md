# Standardize a compound identity audit table

Converts uafR compound-resolution outputs into a stable audit schema for
review, handoff, and reproducible downstream filtering. The audit table
does not invent structures or identifiers; unresolved fields remain
missing.

## Usage

``` r
standardizeCompoundIdentityAudit(x, occurrences = NULL)
```

## Arguments

- x:

  Plant phytochemistry result, \`CompoundResolution\` data frame, or
  compatible identity table.

- occurrences:

  Optional occurrence table used to add query counts when the input
  lacks them.

## Value

Compound identity audit data frame.
