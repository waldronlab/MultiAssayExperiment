context("Subset tests")

example("MultiAssayExperiment")

rowSubsettor <- lapply(rownames(myMultiAssayExperiment)[1:2], function(a) {
  sample(a, 1)
})

test_that("subsettor length is of the same as MultiAssayExperiment", {
  expect_false(length(rowSubsettor) == length(myMultiAssayExperiment))
  expect_error(myMultiAssayExperiment[rowSubsettor, ])
})