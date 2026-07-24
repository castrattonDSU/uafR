# standardifyIt

Post-merge function to quantify compound emission rates relative to
internal or external standard(s). If using internal standard (IS), will
use the contained \`standardify()\` function with user-specified inputs.
If using external standard (ES), an input "matrix" from which standard
curves can be derived.

## Usage

``` r
standardifyIt(
  data_in,
  standard_type = "Internal",
  standard_used = "Tetradecane",
  IS_ng = 190.5,
  IS_uL = 1,
  collect_time = 1,
  sample_amt = 1,
  ES_calibration = NA
)
```

## Arguments

- data_in:

  mzExacto output

- standard_type:

  specifies type of standardization to perform ("Internal" or
  "External")

- standard_used:

  specifies the standard used as an internal standard (name your
  chemical)

- IS_ng:

  specifies the quantified number of molecules in the standard

- IS_uL:

  specifies the amount of standard added to the samples

- collect_time:

  specifies how long samples were collected for (e.g. hour(s), day(s),
  year(s))

- sample_amt:

  specifies how many individuals/items the samples were collected from

- ES_calibration:

  input matrix from which calibrations will be fit for external
  standardization

## Value

Returns the original data standardized relative to the user-specified
internal or external inputs

## Examples

``` r
standardifyIt(standard_exacto, standard_type = "Internal", standard_used = "Octanal",
IS_ng = 1, IS_uL = 1, collect_time = 1)
#>                     Compound          Mass               RT       Best Match
#> 2 Hexanoic acid, ethyl ester 144.115029749 5.37971887422041  99.350118110327
#> 3          Methyl salicylate 152.047344113 8.29568988791984  98.161520884759
#> 4                   Undecane 156.187800766 6.79854308486045 98.6771852019613
#>   Std_soln_00.D Std_soln_07.D Std_soln_00a.D
#> 2    0.69588802     0.3198285       2.345333
#> 3    0.07999628     0.2093027       1.667955
#> 4    0.22751807     0.0000000       1.899181
standardifyIt(standard_exacto, standard_type = "External", ES_calibration = ExternalStandard_data)
#>                     Compound          Mass               RT       Best Match
#> 1                    Octanal 128.120115130  5.4620897525496 99.3245676204613
#> 2 Hexanoic acid, ethyl ester 144.115029749 5.37971887422041  99.350118110327
#> 3          Methyl salicylate 152.047344113 8.29568988791984  98.161520884759
#> 4                   Undecane 156.187800766 6.79854308486045 98.6771852019613
#>   Std_soln_00.D Std_soln_07.D Std_soln_00a.D
#> 1      7.945035      7.936631       7.938918
#> 2      7.942251      7.936124       7.943000
#> 3      7.936617      7.936041       7.940944
#> 4      7.937966            NA       7.941646
```
