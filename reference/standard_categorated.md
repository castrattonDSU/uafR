# Custom Output from \`categorate()\`

A list containing all of the output from
\`categorate(standard_data\$Compound.Name\[standard_data\$Match.Factor
\> 78\])\`. Has all necessary elements for downstream functions,
specifically \`mzExacto()\`.

## Usage

``` r
standard_categorated
```

## Format

\## \`standard_categorated\` A list with 4 data frames.

- Databases:

  Results from Database Searches

- FMCS:

  Results from Atomic Structure Summaries

- FunctionalGroups:

  Results from Structural Similarity Matches with Input Library

- BestChemMatch:

  Library Chemicals that had Strong Structural Matches (Tanimoto
  Similarity \> 0.95)
