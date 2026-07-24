# Enrich plant compounds with existing uafR compound workflows

Enrich plant compounds with existing uafR compound workflows

## Usage

``` r
enrichPlantCompounds(
  plant_compounds,
  chemical_library = NULL,
  detail = c("research", "full", "none"),
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.2,
  enrichment_fun = NULL,
  pubchem_fun = NULL,
  compound_request_fun = NULL,
  batch_size = Inf,
  resume = TRUE,
  progress = interactive(),
  ...
)
```

## Arguments

- plant_compounds:

  A normalized or curated plant-compound table.

- chemical_library:

  Optional library passed to \`categorate()\`.

- detail:

  One of \`"research"\`, \`"full"\`, or \`"none"\`.

- cache:

  Logical.

- cache_dir:

  Cache directory.

- throttle:

  Request throttle.

- enrichment_fun:

  Optional injected enrichment function for tests or cached workflows.

- pubchem_fun:

  Optional replacement for \`pubchemProfile()\` used by the PubChem-only
  fallback. Intended for tests and advanced cached workflows.

- compound_request_fun:

  Optional request function passed to \`pubchemProfile()\` when the
  PubChem-only fallback is used.

- batch_size:

  Maximum number of unique compounds per PubChem-only enrichment batch.
  Use \`Inf\` for one batch.

- resume:

  Logical. If \`TRUE\`, reuse saved PubChem-only batch \`.rds\` files in
  \`cache_dir\`.

- progress:

  Logical. If \`TRUE\`, print simple progress messages during
  PubChem-only enrichment.

- ...:

  Additional arguments passed to \`categorate()\`, \`enrichment_fun\`,
  or the PubChem-only enrichment fallback.

## Value

A list with \`CategorateResult\`, \`CompoundResolution\`,
\`TraitEvidence\`, and \`Provenance\`.
