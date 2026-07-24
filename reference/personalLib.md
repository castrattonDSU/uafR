# Formatting user's personal library

'personalLib' takes a data set (formatted long or wide) and returns
named vectors of the provided data either as separate vectors or a list
of vectors.

## Usage

``` r
personalLib(data, input_Format, output_Format)
```

## Arguments

- data:

  A dataframe or matrix

- input_Format:

  either "long" or "wide"

- output_Format:

  either "vectors" or "list"

## Value

multiple separate character vectors \#' (one per variable type) when
using ouput option "vectors"

a list named "librarylist" containing the multiple separate character
vectors when using output option "list"

## Details

This function allows easy manipulation of an imported data set for
creation into "personal libraries" to use in downstream functions.
Imported data set(s) can be in either long (containing one column for
possible variable types and one column for the values of those variable
types) or wide (containing a column labeled for each variable type with
the values of that variable type in the corresponding rows) format. User
may designate the format of the input data and the desired format for
the output data. Data output options will either be in multiple separate
character vectors (one per variable type) or a list (named
"librarylist") containing the multiple separate character vectors.

## Examples

``` r
personalLib(library_data, "wide", "vectors")
```
