context("drop argument for rows works as it should")

example("MultiAssayExperiment")

test_that("when drop is TRUE, length of Elist is correct", {
  expect_true(length(myMultiAssayExperiment[FALSE, , , drop = TRUE]) == 0L)
  expect_true(length(subset(myMultiAssayExperiment, FALSE, "rows",
                            drop = TRUE)) == 0L)
})
test_that("when drop is FALSE, length is the same", {
  expect_true(length(myMultiAssayExperiment[FALSE, , , drop = FALSE])
              == length(myMultiAssayExperiment))
  expect_true(length(subset(myMultiAssayExperiment, FALSE, "rows",
                            drop = FALSE)) == length(myMultiAssayExperiment))
})

context("drop arugment for columns works as it should")

test_that("when drop is TRUE, length of Elist is correct", {
  expect_true(length(myMultiAssayExperiment[, FALSE, , drop = TRUE]) == 0L)
  expect_true(length(subset(myMultiAssayExperiment, FALSE, "columns",
                            drop = TRUE)) == 0L)
})
test_that("when drop is FALSE, length is the same", {
  expect_true(length(myMultiAssayExperiment[, FALSE, , , drop = FALSE])
              == length(myMultiAssayExperiment))
  expect_true(length(subset(myMultiAssayExperiment, FALSE, "columns",
                            drop = FALSE)) == length(myMultiAssayExperiment))
})