test_that("exactoThese supports bundled Databases schema", {
  reactives = exactoThese(standard_categorated,
                          subsetBy = "Database",
                          subsetArgs = "reactives")
  lotus = exactoThese(standard_categorated,
                      subsetBy = "Database",
                      subsetArgs = "LOTUS")

  expect_true(length(reactives) > 0)
  expect_true("Octanal" %in% reactives)
  expect_true(length(lotus) > 0)
  expect_true("Octanal" %in% lotus)
})

test_that("exactoThese supports split database tables from categorate", {
  categorated = list(
    reactives = data.frame(reactives = c("Reactive Group", "None"),
                           Chemical = c("chem1", "chem2")),
    LOTUS = data.frame(LOTUS = c("L1", "L2"),
                       Chemical = c("chem1", "chem3")),
    KEGG = data.frame(KEGG = "None", Chemical = "chem4"),
    FEMA = data.frame(FEMA = "F1", Chemical = "chem1"),
    FDA_SPL = data.frame(FDA_SPL = "None", Chemical = "chem1")
  )

  expect_equal(exactoThese(categorated, "Database", "reactives"), "chem1")
  expect_equal(exactoThese(categorated, "Database", c("reactives", "LOTUS")),
               "chem1")
  expect_equal(sort(exactoThese(categorated, "Database", "All")),
               c("chem1", "chem3"))
})
