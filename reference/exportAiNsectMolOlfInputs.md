# Export aiNsect molecular-olfaction inputs from wide EO GC-MS tables

\`exportAiNsectMolOlfInputs()\` converts paired wide essential-oil GC-MS
identity and abundance tables into molecular-olfaction input files for
aiNsect. Compound structures are resolved through \`pubchemProfile()\`;
the exporter does not fabricate SMILES, InChIKeys, PubChem CIDs,
formulas, or abundance values.

## Usage

``` r
exportAiNsectMolOlfInputs(
  chem_id_csv,
  chem_quant_csv,
  out_dir,
  cache_dir = NULL,
  profile = "ms",
  throttle = 0.2,
  min_relative_abundance = 0,
  source_label = "uafR_pubchem_eo_gcms",
  uafR_run_id = NULL,
  treatment_aliases = NULL,
  write_pubchem_audit = TRUE,
  profile_fun = pubchemProfile
)
```

## Arguments

- chem_id_csv:

  Path to the wide compound-identity CSV. Columns are treatments and
  rows are ranked compound names.

- chem_quant_csv:

  Path to the matching wide abundance CSV. Columns are treatments and
  rows are abundance values for the same treatment/rank cells.

- out_dir:

  Directory where aiNsect input files will be written.

- cache_dir:

  Directory used by \`pubchemProfile()\` for PubChem response caching.
  If \`NULL\`, \`pubchemProfile()\` uses its default cache directory.

- profile:

  PubChem enrichment profile passed to \`pubchemProfile()\`. Defaults to
  \`"ms"\`.

- throttle:

  Seconds to wait between uncached PubChem requests.

- min_relative_abundance:

  Minimum relative abundance retained per treatment after unresolved
  compounds are removed and duplicate compounds are aggregated. Defaults
  to \`0\`.

- source_label:

  Source label written to \`uafR_compounds.csv\`.

- uafR_run_id:

  Optional run identifier. If \`NULL\`, a timestamped ID is generated.

- treatment_aliases:

  Optional named character vector mapping treatment aliases to canonical
  treatment names. Built-in aliases handle the known EO table cleanup
  issues such as \`"Clary Sagae"\` to \`"Clary Sage"\`.

- write_pubchem_audit:

  Logical. If \`TRUE\`, write raw PubChem identity and property audit
  CSVs when those tables are available.

- profile_fun:

  Advanced/testing hook. Defaults to \`pubchemProfile\`.

## Value

A list containing output paths, exported data frames, unresolved
compounds, and summary metadata. Files written to \`out_dir\` include
\`uafR_compounds.csv\`, \`treatment_compound_abundance.csv\`,
\`uafR_compounds_unresolved.csv\`, and
\`uafR_mololf_export_summary.json\`.

## Examples

``` r
if (FALSE) { # \dontrun{
exportAiNsectMolOlfInputs(
  chem_id_csv = "20240612-EO-gcms-data_all.csv",
  chem_quant_csv = "20240612-EO-gcms-quant_all.csv",
  out_dir = "mololf_export",
  cache_dir = "pubchem_cache",
  profile = "ms",
  throttle = 0.2
)
} # }
```
