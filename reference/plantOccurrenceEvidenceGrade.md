# Grade plant-compound occurrence evidence

\`plantOccurrenceEvidenceGrade()\` turns source evidence fields into a
transparent analysis tier. Grades are conservative and do not claim that
a compound was measured in project samples.

## Usage

``` r
plantOccurrenceEvidenceGrade(x, review_table = NULL)
```

## Arguments

- x:

  Plant phytochemistry result, occurrence table, or enriched
  plant-compound membership table.

- review_table:

  Optional review table created by uafR review workflows.

## Value

Row-level evidence grade table.
