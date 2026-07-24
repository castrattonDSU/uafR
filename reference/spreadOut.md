# spreadOut

Prepares input mass spectrometry data for subsequent functions. Relevant
information for every instance of every chemical across all samples is
placed on the same matrices with no overlaps. Chemicals are also
searched on PubChem and all published names, published m/z histograms,
and exact mass data are downloaded and stored as a large list.

## Usage

``` r
spreadOut(input)
```

## Arguments

- input:

  A .CSV file containing raw mass spectrometry output with
  'Component.RT', 'Component.Area', 'Base.Peak.MZ', 'File.Name',
  'Compound.Name', and 'Match.Factor' as column names/required input

## Value

Returns a list with matrices storing raw area values, chemical
identities, m/z values, match factors, retention times, exact mass data
(if published), and a unique code for each data point (retention time
pasted to exact mass). The list also contains an additional list with
all published chemical names, top m/z peaks, exact mass, and structural
fingerprints in SDF format.

## Details

Automatically prepares .CSV input containing retention time, detected
m/z, component area, compound name, match factor, and sample or file
name for downstream uafR functions.

## Examples

``` r
if (FALSE) { # \dontrun{
spreadOut(standard_data)
} # }
```
