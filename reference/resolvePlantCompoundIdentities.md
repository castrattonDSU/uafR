# Resolve plant compound identities with a fast PubChem pass

Resolves unique compound names from a plant occurrence table or plant
phytochemistry result to PubChem identity fields without pulling the
richer annotation tables used by \`categorate(detail = "research")\`.
This is the preferred first compound-resolution step for large
species-first runs.

## Usage

``` r
resolvePlantCompoundIdentities(
  plant_compounds,
  cache = TRUE,
  cache_dir = NULL,
  throttle = 0.2,
  batch_size = 100,
  resume = TRUE,
  progress = interactive(),
  pubchem_fun = NULL,
  compound_request_fun = NULL,
  lotus_index = NULL,
  source_only = FALSE
)
```

## Arguments

- plant_compounds:

  Plant phytochemistry result, normalized occurrence table, or curated
  plant-compound intake table.

- cache:

  Logical. If \`TRUE\`, PubChem responses and identity batch results can
  be reused.

- cache_dir:

  Cache directory.

- throttle:

  Seconds to wait between uncached PubChem requests.

- batch_size:

  Maximum number of unique compounds per identity batch.

- resume:

  Logical. If \`TRUE\`, reuse saved identity batch \`.rds\` files.

- progress:

  Logical. If \`TRUE\`, print simple progress messages.

- pubchem_fun:

  Optional replacement for \`pubchemProfile()\` used in tests or
  advanced cached workflows.

- compound_request_fun:

  Optional request function passed to \`pubchemProfile()\`.

- lotus_index:

  Optional local LOTUS index as a data frame, flat file, or
  manifest-backed lookup directory. When supplied, source-backed LOTUS
  SMILES, InChIKeys, formulas, and PubChem CIDs are used before PubChem
  name lookup.

- source_only:

  Logical. If \`TRUE\`, resolve only exact identities available from
  \`lotus_index\` and do not make PubChem requests. Unresolved and
  ambiguous names remain explicit for later review or server-side
  resolution.

## Value

A \`CompoundResolution\` data frame. The underlying identity-only
PubChem tables are attached as the \`"CategorateResult"\` attribute.
Compound keys preserve chemically meaningful Greek-letter and plus/minus
prefixes so isomers such as alpha-pinene and beta-pinene are not
collapsed together.
