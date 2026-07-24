test_that("external standardization inverts the selected calibration model", {
  old_scipen = getOption("scipen")
  standardized = standardifyIt(standard_exacto,
                               standard_type = "External",
                               ES_calibration = ExternalStandard_data)

  ES_Abundance = ExternalStandard_data$Component.Area
  ES_ng = ExternalStandard_data$Quantity
  models = list(log = lm(ES_Abundance ~ log(ES_ng)),
                exponent = lm(ES_Abundance ~ exp(ES_ng)),
                linear = lm(ES_Abundance ~ ES_ng))
  r_squared = vapply(models, function(model) {
    summary(model)$adj.r.squared
  }, numeric(1))
  best_model_name = names(r_squared)[which.max(r_squared)]
  best_model = models[[best_model_name]]
  coefficients = as.numeric(best_model$coefficients)
  current_area = as.numeric(standard_exacto[1, 5])

  expected = switch(best_model_name,
                    log = exp((current_area - coefficients[1]) / coefficients[2]),
                    exponent = log((current_area - coefficients[1]) / coefficients[2]),
                    linear = (current_area - coefficients[1]) / coefficients[2])

  expect_equal(as.numeric(standardized[1, 5]), expected)
  expect_equal(getOption("scipen"), old_scipen)
})

test_that("standardifyIt rejects ambiguous standardization inputs", {
  expect_error(standardifyIt(standard_exacto, standard_type = "Bogus"),
               "standard_type")
  expect_error(standardifyIt(standard_exacto,
                             standard_type = "Internal",
                             standard_used = "No such chemical"),
               "standard_used")
})
