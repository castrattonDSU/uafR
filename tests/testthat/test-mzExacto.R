
run_live_uafr_tests = function() {
 identical(Sys.getenv("UAFR_RUN_LIVE_TESTS"), "true")
}


test_that("output is always correct size", {
  skip_if_not(run_live_uafr_tests(),
              "Live PubChem/NCI integration test")
  search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("ethyl hexanoate", "undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))

  search_chems = c("undecane")
  expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(search_chems))
})

test_that("duplicates do not matter", {
 skip_if_not(run_live_uafr_tests(),
             "Live PubChem/NCI integration test")
 search_chems = c("ethyl hexanoate", "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "ethyl hexanoate", "methyl salicylate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane", "methyl salicylate", "octanal")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane",
                  "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_equal(nrow(mzExacto(standard_spread, search_chems)), length(unique(search_chems)))
})

test_that("missing query chemicals are bad",{
 skip_if_not(run_live_uafr_tests(),
             "Live PubChem/NCI integration test")
 search_chems = c("", "ethyl hexanoate", "methyl salicylate", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("", "", "methyl salicylate", "octanal", "undecane")
 expect_error(mzExacto(standard_spread, search_chems))

 search_chems = c("ethyl hexanoate", "methyl salicylate", "octanal", "undecane", "", "")
 expect_error(mzExacto(standard_spread, search_chems))
})
