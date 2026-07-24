# mzExacto

Uses the output from \`spreadOut()\` and a list of chemicals to extract
data for. Samples that contain a chemical have all identified area(s)
aggregated. The most likely identifications are determined by
prioritizing matches of exact chemical names, followed by m/z overlaps
within precise retention time windows that are determined by molecular
masses published on PubChem.

## Usage

``` r
mzExacto(data_in, chemicals, decontaminate = T)
```

## Arguments

- data_in:

  the list output from \`spreadOut()\`

- chemicals:

  A vector containing chemical names in IUPAC notation

- decontaminate:

  A logical set to remove uncertain chemicals or junk from output

## Value

A data frame containing chemical names(\`chemicals\`), their optimal
retention time, exact mass, best identified match factor, and aggregated
component area across every sample it was identified in.

## Details

Communicates with PubChem to collect information on every search
chemical. Uses this information to search the raw \`spreadOut()\` data
for samples where a chemical exists.

## Examples

``` r
if (FALSE) { # \dontrun{
query_chemicals = c("Ethyl hexanoate","Methyl salicylate","Octanal","Undecane")
mzExacto(standard_spread, query_chemicals)
} # }
```
