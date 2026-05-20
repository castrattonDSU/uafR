
run_live_uafr_tests = function() {
  identical(Sys.getenv("UAFR_RUN_LIVE_TESTS"), "true")
}

query_chemicals = c("Linalool", "alpha-Pinene", "Aspirin", "Caffeine", "Limonene")

test_that("categorate requires compounds and a chemical library", {
  expect_error(categorate(chemical_library = library_data),
               "`compounds` is required")
  expect_error(categorate(query_chemicals),
               "`chemical_library` is required")
})

test_that("output has the correct number of objects",{
  skip_if_not(run_live_uafr_tests(),
              "Live PubChem/NCI integration test")
  query_categorated = suppressWarnings(categorate(query_chemicals, library_data, input_format = "wide"))
  expect_equal(length(query_categorated), 8)
  expect_true(all(c("reactives", "LOTUS", "KEGG", "FEMA", "FDA_SPL",
                    "FMCS", "FunctionalGroups", "BestChemMatch") %in%
                    names(query_categorated)))
})

test_that("only items in output are for query chemicals",{
 skip_if_not(run_live_uafr_tests(),
             "Live PubChem/NCI integration test")
 query_categorated = suppressWarnings(categorate(query_chemicals, library_data, input_format = "wide"))
 expect_true(all(unique(query_categorated$reactives$Chemical) %in% query_chemicals))
 expect_equal(length(unique(query_categorated$FMCS$Chemical)), length(query_chemicals))
 expect_equal(length(unique(query_categorated$FunctionalGroups$Chemical)), length(query_chemicals))
 expect_equal(length(unique(query_categorated$BestChemMatch$Chemical)), length(query_chemicals))
})
