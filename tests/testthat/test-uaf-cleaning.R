test_that("CID extraction handles webchem service fallback shapes", {
  expect_equal(.uaf_extract_cid(data.frame(query = "aspirin", cid = 2244)),
               "2244")
  expect_true(is.na(.uaf_extract_cid(NA)))
  expect_true(is.na(.uaf_extract_cid(data.frame())))
  expect_true(is.na(.uaf_extract_cid("Limit Met")))
})
