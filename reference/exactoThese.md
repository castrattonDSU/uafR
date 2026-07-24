# exactoThese

Takes categorated output as input and makes chemical subsets by the
information in it.

## Usage

``` r
exactoThese(
  categoratedInput,
  subsetBy = "Database",
  subsetArgs = "All",
  subsetArgs2 = NA,
  subset_input = NA
)
```

## Arguments

- categoratedInput:

  the lists output from \`categorate()\`

- subsetBy:

  specifies the list to subset by (Database, FMCS, or Library)

- subsetArgs:

  additional arguments to specify which values to subset by

- subsetArgs2:

  additional arguments to specify which values to subset by

- subset_input:

  used when subsetting by FMCS information (e.g. molecular weight)

## Value

A vector of chemical names that meet each user-specification. Used as
input for \`mzExacto()\`.

## Details

Provides a set of search chemicals for \`mzExacto()\`. User gets to
select chemicals based on information generated from \`categorate()\`.
Database subsetting accepts both historical \`Databases\` output and
current split database tables from \`categorate()\`.

## Examples

``` r
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "All")
#>  [1] "Boron trifluoride"                   
#>  [2] "Methylene chloride"                  
#>  [3] "Borane carbonyl"                     
#>  [4] "Octanal"                             
#>  [5] "Undecane"                            
#>  [6] "Methyl salicylate"                   
#>  [7] "Acetamide, 2-fluoro-"                
#>  [8] "Hexanoic acid, ethyl ester"          
#>  [9] "2-Hexen-1-ol, (E)-"                  
#> [10] "2-Octen-1-ol, (E)-"                  
#> [11] "Benzene, 1,3-bis(1,1-dimethylethyl)-"
#> [12] "1H-Tetrazol-5-amine"                 
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "reactives")
#> [1] "Boron trifluoride"    "Methylene chloride"   "Borane carbonyl"     
#> [4] "Octanal"              "Undecane"             "Methyl salicylate"   
#> [7] "Acetamide, 2-fluoro-"
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "LOTUS")
#> [1] "Methylene chloride"                  
#> [2] "Borane carbonyl"                     
#> [3] "Hexanoic acid, ethyl ester"          
#> [4] "Octanal"                             
#> [5] "2-Hexen-1-ol, (E)-"                  
#> [6] "2-Octen-1-ol, (E)-"                  
#> [7] "Undecane"                            
#> [8] "Methyl salicylate"                   
#> [9] "Benzene, 1,3-bis(1,1-dimethylethyl)-"
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "KEGG")
#> [1] "Borane carbonyl"      "Octanal"              "Methyl salicylate"   
#> [4] "Acetamide, 2-fluoro-"
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FEMA")
#> [1] "Borane carbonyl"            "Hexanoic acid, ethyl ester"
#> [3] "2-Hexen-1-ol, (E)-"         "Methyl salicylate"         
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = "FDA_SPL")
#>  [1] "Boron trifluoride"          "Methylene chloride"        
#>  [3] "Borane carbonyl"            "Hexanoic acid, ethyl ester"
#>  [5] "Octanal"                    "2-Hexen-1-ol, (E)-"        
#>  [7] "2-Octen-1-ol, (E)-"         "Undecane"                  
#>  [9] "1H-Tetrazol-5-amine"        "Methyl salicylate"         
#> [11] "Acetamide, 2-fluoro-"      
exactoThese(standard_categorated, subsetBy = "Database", subsetArgs = c("reactives", "FEMA"))
#> [1] "Borane carbonyl"   "Methyl salicylate"
exactoThese(standard_categorated, subsetBy = "FMCS", subsetArgs = "MW",
subsetArgs2 = "Between", subset_input = c(125, 200))
#> [1] "hexanoic acid ethyl ester"          "octanal"                           
#> [3] "(E)-2-octen-1-ol"                   "undecane"                          
#> [5] "2,2,3,4-tetramethylpentane"         "2-hydroxybenzoic acid methyl ester"
#> [7] "1,3-ditert-butylbenzene"           
```
