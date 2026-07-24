test_that("CID extraction handles webchem service fallback shapes", {
  expect_equal(.uaf_extract_cid(data.frame(query = "aspirin", cid = 2244)),
               "2244")
  expect_true(is.na(.uaf_extract_cid(NA)))
  expect_true(is.na(.uaf_extract_cid(data.frame())))
  expect_true(is.na(.uaf_extract_cid("Limit Met")))
})

test_that("yes/no helper preserves vectorized row-level decisions", {
  expect_equal(.uaf_yes_no(TRUE), "Yes")
  expect_equal(.uaf_yes_no(FALSE), "No")
  expect_equal(.uaf_yes_no(c(TRUE, FALSE, NA)), c("Yes", "No", "No"))
})
