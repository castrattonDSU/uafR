# uafR

## Automated GC/LC-MS data processing

------------------------------------------------------------------------

Nothing in life is to be feared; it is only to be understood.  
–Marie Curie

------------------------------------------------------------------------

![Graphical Abstract](GraphicalAbstract.jpg)

Graphical Abstract

### Hello \[Chemical\] World!

       Chemistry plays an active role in every aspect of human
existence. Whether it be the digested molecules of our food, the solid
matrices of plastic and metal that comprise our technology, or the
macromolecules that build and run our cells, the compositions are all
chemical. To understand this aspect of our existence, there are advanced
instruments and techniques that identify chemicals of any system to
name. These precise instruments also allow the number of molecules for
each individual chemical to be quantified across analyzed samples. While
the output from these machines is immediately available, preparing the
raw output for interpretive statistics can require hours/days of trained
labor per sample and months per experiment.   
  
       To remove the bottleneck that exists between the acquisition of
raw mass spectrometry output and the interpretation of chemicals across
experimental treatments, we have developed advanced algorithms that
automate the entire process. Ours is the first GC/LC-MS utility that
accesses published information for every tentative compound to
intelligently select portions of each sample that describe
user-specified query chemicals. This information makes uafR the most
accurate and advanced post GC/LC-MS processing application to date.

### New Standard, Mass Spectrometry Workflows

#### Input Data Structure

       The original workflow for uafR was developed using Agilent
instruments and software. The recommended software for generating the
necessary data in the default format (i.e. with correct column names) is
[Unknowns
Analysis](https://www.agilent.com/cs/library/usermanuals/public/G3335-90187_Unknowns_Analysis_Familiarization-en.pdf).
That said, any software or utility that generates the necessary
information can be used with simple modifications (e.g. changing the
column names).

| Component.RT | Base.Peak.MZ | Component.Area | Compound.Name | Match.Factor | File.Name |
|:--:|:--:|:--:|:---|:--:|:--:|
| 8.229034 | 84.00 | 906.4701 | Pipradrol | 62.62271 | Std_soln_07 |
| 8.286703 | 120.00 | 209705.1878 | Methyl salicylate | 98.16152 | Std_soln_00a |
| 8.296408 | 119.99 | 30332.9022 | Methyl salicylate | 95.79911 | Std_soln_00 |
| 8.303958 | 120.00 | 6476.4785 | Methyl salicylate | 86.29569 | Std_soln_07 |
| 8.348031 | 105.00 | 420.8119 | 3-Hexen-1-ol, benzoate, (Z)- | 68.78156 | Std_soln_00 |
| **…** | **…** | **…** | **…** | **…** | **…** |

#### Spread It Out

       The first step in the process is to convert the raw input to a
format that downstream functions can work with.
[`spreadOut()`](https://castrattonDSU.github.io/uafR/reference/spreadOut.md)
prepares the read in .CSV for intelligent ***sorting*** (using retention
times and published masses) then ***aggregation*** (using all published
names and top m/z peaks) of sample portions that describe a chemical. A
list containing all necessary information for the next function,
[`mzExacto()`](https://castrattonDSU.github.io/uafR/reference/mzExacto.md),
is returned.  
       Contents of the list include matrices (here focused on methyl
salicylate) that store:

1.  **chemical names**

|           Sample_1           |     Sample_2      |     Sample_3      |
|:----------------------------:|:-----------------:|:-----------------:|
|            \<NA\>            |     Pipradrol     |      \<NA\>       |
|            \<NA\>            |      \<NA\>       | Methyl salicylate |
|      Methyl salicylate       |      \<NA\>       |      \<NA\>       |
|            \<NA\>            | Methyl salicylate |      \<NA\>       |
| 3-Hexen-1-ol, benzoate, (Z)- |      \<NA\>       |      \<NA\>       |
|            **…**             |       **…**       |       **…**       |

2.  **retention times**

|  Sample_1   |  Sample_2   |  Sample_3   |
|:-----------:|:-----------:|:-----------:|
|   \<NA\>    | 8.229033559 |   \<NA\>    |
|   \<NA\>    |   \<NA\>    | 8.286703432 |
| 8.296408204 |   \<NA\>    |   \<NA\>    |
|   \<NA\>    | 8.303958027 |   \<NA\>    |
| 8.348031108 |   \<NA\>    |   \<NA\>    |
|    **…**    |    **…**    |    **…**    |

3.  **match factors**

|  Sample_1   |  Sample_2   |  Sample_3   |
|:-----------:|:-----------:|:-----------:|
|   \<NA\>    | 62.62271472 |   \<NA\>    |
|   \<NA\>    |   \<NA\>    | 98.16152088 |
| 95.79911297 |   \<NA\>    |   \<NA\>    |
|   \<NA\>    | 86.29569222 |   \<NA\>    |
| 68.78156469 |   \<NA\>    |   \<NA\>    |
|    **…**    |    **…**    |    **…**    |

4.  **captured M/Z value**

| Sample_1 | Sample_2 | Sample_3 |
|:--------:|:--------:|:--------:|
|  \<NA\>  |    84    |  \<NA\>  |
|  \<NA\>  |  \<NA\>  |   120    |
|  119.99  |  \<NA\>  |  \<NA\>  |
|  \<NA\>  |   120    |  \<NA\>  |
|   105    |  \<NA\>  |  \<NA\>  |
|  **…**   |  **…**   |  **…**   |

5.  **exact mass data (if published)**

|   Sample_1    |   Sample_2    |   Sample_3    |
|:-------------:|:-------------:|:-------------:|
|    \<NA\>     | 267.162314293 |    \<NA\>     |
|    \<NA\>     |    \<NA\>     | 152.047344113 |
| 152.047344113 |    \<NA\>     |    \<NA\>     |
|    \<NA\>     | 152.047344113 |    \<NA\>     |
| 204.115029749 |    \<NA\>     |    \<NA\>     |
|     **…**     |     **…**     |     **…**     |

6.  **raw area values**

|  Sample_1   |  Sample_2   |  Sample_3   |
|:-----------:|:-----------:|:-----------:|
|   \<NA\>    | 906.4700739 |   \<NA\>    |
|   \<NA\>    |   \<NA\>    | 209705.1878 |
| 30332.90221 |   \<NA\>    |   \<NA\>    |
|   \<NA\>    | 6476.478451 |   \<NA\>    |
| 420.8119135 |   \<NA\>    |   \<NA\>    |
|    **…**    |    **…**    |    **…**    |

7.  **a unique code for each input data point (retention time pasted to
    exact mass)**

| Sample_1 | Sample_2 | Sample_3 |
|:--:|:--:|:--:|
| \<NA\> | 8.229033559 \| 267.162314293 | \<NA\> |
| \<NA\> | \<NA\> | 8.286703432 \| 152.047344113 |
| 8.296408204 \| 152.047344113 | \<NA\> | \<NA\> |
| \<NA\> | 8.303958027 \| 152.047344113 | \<NA\> |
| 8.348031108 \| 204.115029749 | \<NA\> | \<NA\> |
| **…** | **…** | **…** |

8.  **and a nested list with**:

- *all published chemical names* (only first 5 are shown)

&nbsp;

    #> [1] "methyl salicylate"        "Methyl 2-hydroxybenzoate"
    #> [3] "119-36-8"                 "Wintergreen oil"         
    #> [5] "Gaultheria oil"

- *top m/z peaks*

&nbsp;

    #> [1] "120" "92"  "152" "121" "65"

- *exact mass*

&nbsp;

    #> [1] 152.0473

- *and likely retention times for the query chemicals.*

&nbsp;

    #> [1] 8.286703

#### Extract Your Chemicals

       The output from
[`spreadOut()`](https://castrattonDSU.github.io/uafR/reference/spreadOut.md)
is like a searchable chemical database where each entry has every
published, uniquely identifying feature assigned to it.
[`mzExacto()`](https://castrattonDSU.github.io/uafR/reference/mzExacto.md)
collects the same information for a set of query chemicals and uses it
to precisely search the advanced dictionary for samples that have those
chemicals.   
       In many cases, users will already know what they are looking for.
In others, they won’t.

------------------------------------------------------------------------

##### When Chemicals are Known 

       While there are multiple ways to create a list of input chemicals
\[see
[`personalLib()`](https://castrattonDSU.github.io/uafR/reference/personalLib.md)\],
a simple method for smaller searches is to just type quotes around the
search names in a list:

**`query_chemicals = c("Ethyl hexanoate", "Methyl salicylate", "Octanal", "Undecane")`**

[`mzExacto()`](https://castrattonDSU.github.io/uafR/reference/mzExacto.md)
takes the output from spreadOut() \[`standard_spread`\] and this list of
`query_chemicals`:

**`mzExacto(standard_spread, query_chemicals)`**

returning a single dataframe with all of the necessary information for
downstream functions and, ultimately, interpretation.

| Compound | Mass | RT | Best Match | Std_soln_00 | Std_soln_07 | Std_soln_00a |
|----|----|----|----|----|----|----|
| Octanal | 128.120115130 | 5.462089753 | 99.32456762 | 379178.88653 | 30943.11385 | 125725.8982 |
| Ethyl hexanoate | 144.115029749 | 5.379718874 | 99.35011811 | 263866.0427 | 9896.488149 | 294869.1357 |
| Methyl salicylate | 152.047344113 | 8.295689887 | 98.16152088 | 30332.90221 | 6476.478451 | 209705.1878 |
| Undecane | 156.187800766 | 6.129191467 | 98.6771852 | 86270.05019 | 243.9123731 | 238776.2287 |

------------------------------------------------------------------------

##### When Types/Classes of Chemicals are Known 

       [`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
is an overpowered function that accesses a broad array of categorical
data for searched chemicals. Here we present a single application from
the output of
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
that could help in GC/LC-MS analyses where classes/types of chemical
groupings are known, but not specific compounds. For a detailed overview
of additional chemistry workflows (e.g. meta-analyses) that
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
could catalyze, we recommend the companion manuscript (linked when
published).  
       A required input for running
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
is a library to perform structural matches against. This library is a
.CSV file that can contain as many sets of chemicals as the user’s
hardware can handle. For the following example, we have restricted our
search to 4 sets of chemicals that we know are structurally similar to
our query chemicals from the previous section (Types B, C, D, and E
below) and one set that should not have any matches (Type A).

| Type.A | Type.B | Type.C | Type.D | Type.E |
|:---|:---|:---|:---|:---|
| 2-Aminothiazole | o-Cresol | Nonane | Butyl methacrylate | Octane |
| 3,4,5,6-Tetrachlorocyclohexene | Salicylic Acid | Dodecane | Isobutyl hexanoate | Octanoic acid |
| N-methyl-1,3,5-triazin-2-amine | Guaiacol | Tridecane | Ethyl heptanoate | 1-Octanol |
| 2-Methyloctahydro-2-azacyclopropa\[cd\]pentalene | Aspirin |  | Methyl hexanoate | Hexadecanal |
|  | Salicyl alcohol |  | Methyl heptanoate | Decanal |
|  | 4-Methylsalicylic acid |  | Dihexyl adipate | Undecanal |
|  |  |  | 2-Heptanone | Hexyl acetate |

        To perform the chemical structure matches and summarize atomic
features, uafR taps into an amazing set of cheminformatics packages –
[ChemmineR](https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html),
[fmcsR](https://bioconductor.org/packages/release/bioc/html/fmcsR.html),
[webchem](https://cran.r-project.org/web/packages/webchem/index.html).
The library tests return the following data frames:

| Type.A | Type.B | Type.C | Type.D | Type.E |          Chemical |
|:------:|:------:|:------:|:------:|:------:|------------------:|
|   No   |   No   |   No   |  Yes   |   ~    |   ethyl hexanoate |
|   No   |  Yes   |   No   |   No   |   No   | methyl salicylate |
|   No   |   No   |   ~    |   ~    |  Yes   |           octanal |
|   No   |   No   |  Yes   |   ~    |  Yes   |          undecane |

       Where “No” means none of the chemicals had a structural match,
“~” refers to at least 1 match between 0.85 and 0.95, and “Yes” means
there was at least 1 match exceeding 0.95.

| Type.A | Type.B | Type.C | Type.D | Type.E |          Chemical |
|:------:|:------:|:------:|:------:|:------:|------------------:|
|   No   |   No   |   No   |  CMP2  |   No   |   ethyl hexanoate |
|   No   |  CMP1  |   No   |   No   |   No   | methyl salicylate |
|   No   |   No   |   No   |   No   |  CMP1  |           octanal |
|   No   |   No   |  CMP1  |   No   |  CMP1  |          undecane |

       Where the first compound in a set that had a match exceeding 0.95
is shown. The number following “CMP” refers tells the user which
compound was a match (i.e. 1 refers to the topmost chemical in the
group), so ethyl hexanoate was more structurally similar to isobutyl
hexanoate than butyl methacrylate. Makes sense!  
       As can be seen, our library did a great job pulling out the
chemicals of interest from the previous example. However, some studies
may have even less direction to go off of. In these cases, the
atomic/functional group summary provided by
[fmcsR](https://bioconductor.org/packages/release/bioc/html/fmcsR.html)
can also help navigate:

| Chemical          | Groups | GroupCounts | Atom | AtomCounts | NCharges |
|:------------------|:------:|:-----------:|:----:|:----------:|:--------:|
| ethyl hexanoate   | RCOOH  |      0      |  H   |     16     |    0     |
| ethyl hexanoate   | RCOOR  |      1      |  O   |     2      |    0     |
| ethyl hexanoate   |  ROR   |      0      |  C   |     8      |    0     |
| methyl salicylate |  ROH   |      1      |  H   |     8      |    0     |
| methyl salicylate | RCOOR  |      1      |  O   |     3      |    0     |
| methyl salicylate |  RCOR  |      0      |  C   |     8      |    0     |
| octanal           |  ROH   |      0      |  H   |     16     |    0     |
| octanal           |  RCHO  |      1      |  O   |     1      |    0     |
| octanal           |  RCOR  |      0      |  C   |     8      |    0     |
| undecane          |  RCHO  |      0      |  H   |     24     |    0     |
| undecane          |  RCOR  |      0      |  C   |     11     |    0     |

       This is only a subset of the output from this utility/function
meant to emphasize some of the more useful results. The actual output
includes columns for the molecular weight, molecular formula,
presence/absence of rings (cyclical carbon groups), and additional
common functional groups including those with phosporous or nitrogen.   
       Getting back to the GC/LC-MS workflow, subsetting the
`query_chemicals` with these outputs is very easily achieved using
[`exactoThese()`](https://castrattonDSU.github.io/uafR/reference/exactoThese.md):

**`query_chemicals = exactoThese(input_categorated, subsetBy = "FMCS",`**  
**`subsetArgs = "MW", subsetArgs2 = "Between", subset_input = c(50,115))`**

       With these arguments,
[`exactoThese()`](https://castrattonDSU.github.io/uafR/reference/exactoThese.md)
returns every input chemical with a molecular weight between 50 and 115
g/mol.

------------------------------------------------------------------------

##### Unknown Exploration

       As referenced,
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
casts a broad net for categorical information on query chemicals. This
internet…net uses regular expressions to pull published data from a
variety of databases. Captured information can include:

1. Reactive Groups from [PubChem](https://pubchem.ncbi.nlm.nih.gov/),

    #>                         reactives_df          Chemical
    #> 41                           Octanal           Octanal
    #> 42                         Aldehydes           Octanal
    #> 43                   N-OCTYLALDEHYDE           Octanal
    #> 48                          Undecane          Undecane
    #> 49 Hydrocarbons, Aliphatic Saturated          Undecane
    #> 60                 Methyl Salicylate Methyl salicylate
    #> 62               Phenols and Cresols Methyl salicylate
    #> 63             PubChem Internal Link Methyl salicylate
    #> 64                 METHYL SALICYLATE Methyl salicylate

2. natural products occurrences from
[LOTUS](https://lotus.naturalproducts.net/),

    #>    LOTUS_df          Chemical
    #> 41  Q416673           Octanal
    #> 42  biochem           Octanal
    #> 48  Q150731          Undecane
    #> 49  biochem          Undecane
    #> 60  Q407669 Methyl salicylate
    #> 61  biochem Methyl salicylate

3. bioactivites and risk categories from the Kyoto Encyclopedia of Genes
and Genomes ([KEGG](https://www.genome.jp/kegg/)),

    #>                                      KEGG_df          Chemical
    #> 41                               KEGG: Lipid           Octanal
    #> 48                                      None          Undecane
    #> 60                                KEGG: Drug Methyl salicylate
    #> 61                                KEGG: JP15 Methyl salicylate
    #> 62 KEGG: Risk Category of Japanese OTC Drugs Methyl salicylate
    #> 63                           KEGG: OTC drugs Methyl salicylate
    #> 64                        KEGG: Animal Drugs Methyl salicylate
    #> 65                         KEGG: Drug Groups Methyl salicylate
    #> 66                        KEGG: Drug Classes Methyl salicylate
    #> 67                                      <NA> Methyl salicylate

4. flavors, odors, etc. from the Flavor and Extract Manufacturers
Association ([FEMA](https://www.femaflavor.org/)),

    #>        FEMA_df          Chemical
    #> 41        None           Octanal
    #> 48        None          Undecane
    #> 60      Almond Methyl salicylate
    #> 61     Caramel Methyl salicylate
    #> 62  Peppermint Methyl salicylate
    #> 63       Sharp Methyl salicylate

5. and whether it exists in the Food and Drug Administration’s SPL data
base ([FDA/SPL](https://www.fda.gov/)).

    #>           FDA_SPL_df          Chemical
    #> 41    CAPRYLALDEHYDE           Octanal
    #> 48          UNDECANE          Undecane
    #> 60 METHYL SALICYLATE Methyl salicylate

       Again, we can easily subset our `query_chemicals()` with this
information using
[`exactoThese()`](https://castrattonDSU.github.io/uafR/reference/exactoThese.md):

**`query_chemicals = exactoThese(chems_categorated, subsetBy = "Database", subsetArgs = c("LOTUS", "FEMA"))`**

       Here,
[`exactoThese()`](https://castrattonDSU.github.io/uafR/reference/exactoThese.md)
returns every input chemical for which information could be found on
both LOTUS and FEMA. While the following output shows the chemicals from
the first search (known chemicals), this example is based on the
assumption that we do not know what we will find across every sample.
The compounds are the focus for this portion of the example only for
clarity and simplicity. In this alternate, unknown context, a useful
approach for narrowing the search chemicals for
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
and/or
[`mzExacto()`](https://castrattonDSU.github.io/uafR/reference/mzExacto.md)
is to first subset by match factor:

**`query_chems = standard_dat$Compound.Name[standard_dat$Match.Factor >= 65]`**

       At this match factor, the example input data structure would
change to:

| Component.RT | Base.Peak.MZ | Component.Area | Compound.Name | Match.Factor | Sample.Name |
|:--:|:--:|:--:|:---|:--:|:--:|
| \<NA\> | \<NA\\ | \<NA\> | \<NA\> | \<NA\> | \<NA\> |
| 8.286703 | 120.00 | 209705.1878 | Methyl salicylate | 98.16152 | Std_soln_00a |
| 8.296408 | 119.99 | 30332.9022 | Methyl salicylate | 95.79911 | Std_soln_00 |
| 8.303958 | 120.00 | 6476.4785 | Methyl salicylate | 86.29569 | Std_soln_07 |
| 8.348031 | 105.00 | 420.8119 | 3-Hexen-1-ol, benzoate, (Z)- | 68.78156 | Std_soln_00 |
| **…** | **…** | **…** | **…** | **…** | **…** |

       Which still leaves some “junk” that our analysis would probably
be better without. To remedy this, we could - 1)
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
all of the chemicals (`Compound.Name`) at this `Match.Factor` then
subset by some feature(s); or, 2) continue to adjust the `Match.Factor`
until the data make more chemical sense.

**`query_chems = standard_dat$Compound.Name[standard_dat$Match.Factor >= 80]`**

       This match factor yields:

| Component.RT | Base.Peak.MZ | Component.Area | Compound.Name | Match.Factor | Sample.Name |
|:--:|:--:|:--:|:---|:--:|:--:|
| \<NA\> | \<NA\\ | \<NA\> | \<NA\> | \<NA\> | \<NA\> |
| 8.286703 | 120.00 | 209705.1878 | Methyl salicylate | 98.16152 | Std_soln_00a |
| 8.296408 | 119.99 | 30332.9022 | Methyl salicylate | 95.79911 | Std_soln_00 |
| 8.303958 | 120.00 | 6476.4785 | Methyl salicylate | 86.29569 | Std_soln_07 |
| \<NA\> | \<NA\\ | \<NA\> | \<NA\> | \<NA\> | \<NA\> |
| **…** | **…** | **…** | **…** | **…** | **…** |

       While it may seem as though all is well from this fraction of the
data, it is important to remember that it is only a peek at what is lost
or gained from the `Match.Factor` adjustments. To emphasize this point,
consider what is lost by this adjustment:

**`query_chems = standard_dat$Compound.Name[standard_dat$Match.Factor >= 90]`**

       The high level of stochasticity behind every data point in a mass
spectrometry analysis is another reason previous algorithms fail when
assigning area values across samples. ***Hopefully*** more simply put –
with chemicals, things don’t always go exactly the same. While this has
historically driven manual over data-driven workflows, it could also be
argued that it drives subjectivity into chemical analysis. Modern
programming languages allow even complex workflows to be automated. By
accessing published information we are able to mirror the optimal manual
workflow for chemicals that were either misread in a sample or buried by
similarities. This allows razor-sharp precision when excising from
“gray” regions of the data.   
       In this example, the known chemicals were found simply by
sub-setting with `Match.Factor`:

**`query_chems = standard_dat$Compound.Name[standard_dat$Match.Factor > 89]`**

**`mzExacto(standard_spread, query_chems)`**

| Compound | Mass | RT | Best Match | Std_soln_00 | Std_soln_07 | Std_soln_00a |
|:---|---:|---:|---:|---:|---:|---:|
| Octanal | 128.120115130 | 5.462089753 | 99.32456762 | 379178.88653 | 30943.11385 | 125725.8982 |
| Ethyl hexanoate | 144.115029749 | 5.379718874 | 99.35011811 | 263866.0427 | 9896.488149 | 294869.1357 |
| Methyl salicylate | 152.047344113 | 8.295689887 | 98.16152088 | 30332.90221 | 6476.478451 | 209705.1878 |
| Undecane | 156.187800766 | 6.129191467 | 98.6771852 | 86270.05019 | 243.9123731 | 238776.2287 |

       But, the combined `Match.Factor` and
[`categorate()`](https://castrattonDSU.github.io/uafR/reference/categorate.md)
approach can churn through a large amount of complex chemical data
faster and with more accuracy than any manual protocol for unknown
compound selections.

------------------------------------------------------------------------
