context("Subset tests")

example("MultiAssayExperiment")

colList <- colnames(myMultiAssayExperiment)
colList[[2L]] <- character(0L)

rowList <- rownames(myMultiAssayExperiment)
rowList[[3L]] <- character(0L)

rowSubsettor <- lapply(rownames(myMultiAssayExperiment)[1:2], function(a) {
  sample(a, 1)
})
subsettor2 <- rownames(myMultiAssayExperiment)

test_that("subsettor length is of the same as MultiAssayExperiment", {
  expect_false(length(rowSubsettor) == length(myMultiAssayExperiment))
  expect_error(myMultiAssayExperiment[rowSubsettor, ])
  expect_equal(myMultiAssayExperiment[subsettor2, ], myMultiAssayExperiment)
})

test_that("drop argument works", {
  expect_equal(length(
      experiments(myMultiAssayExperiment[, colList, drop = TRUE])), 2L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[, colList, drop = FALSE])), 3L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[rowList, drop = TRUE])), 2L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[rowList, drop = FALSE])), 3L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[FALSE, drop = TRUE])), 0L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[FALSE, drop = FALSE])), 3L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[,FALSE, drop = TRUE])), 0L)
  expect_equal(length(
      experiments(myMultiAssayExperiment[,FALSE, drop = FALSE])), 3L)
})
