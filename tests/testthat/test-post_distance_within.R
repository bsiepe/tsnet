test_that("post_distance_within calculates distances correctly for empirical models", {
  data("ts_data")
  data1 <- subset(ts_data, id == "ID1")
  data1 <- subset(data1, select = -c(id))
  fitobj <- BGGM::var_estimate(data1)
  comp <- "frob"
  draws <- 1000

  result <- post_distance_within(fitobj, comp, pred = FALSE, draws = draws)

  # Check if the result has the correct number of rows
  expect_equal(nrow(result), draws)

  # Check if the calculated distances are numeric
  expect_true(all(is.na(result$beta) | is.numeric(result$beta)))
  expect_true(all(is.na(result$pcor) | is.numeric(result$pcor)))

})
