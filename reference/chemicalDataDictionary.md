# Describe categorate output schemas

`chemicalDataDictionary()` returns the expected table and column
contracts for the main `categorate(detail = "research")` and
`categorate(detail = "full")` outputs. It is intended for users who need
to understand, export, validate, or join enriched chemical data without
guessing what each table contains.

## Usage

``` r
chemicalDataDictionary(tables = NULL)
```

## Arguments

- tables:

  Optional character vector of table names to return. If `NULL`, all
  documented enriched output tables are returned.

## Value

A data frame with table names, column names, expected types, required
flags, semantic roles, allowed values, and descriptions.

## Examples

``` r
if (FALSE) { # \dontrun{
dictionary = chemicalDataDictionary()
chemicalDataDictionary("ChemicalTraitReport")
} # }
```
