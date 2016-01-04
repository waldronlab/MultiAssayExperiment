context("Endomorphic tests of builder functions")

test_that("builders return empty", {
  expect_true(isEmpty(MultiAssayExperiment()))
  expect_true(isEmpty(RangedRaggedAssay()))
  expect_true(isEmpty(Elist()))
})

test_that("builders return approrpriate functions", {
  expect_equal(class(MultiAssayExperiment()), "MultiAssayExperiment")
  expect_equal(class(RangedRaggedAssay()), "RangedRaggedAssay")
  expect_equal(class(Elist()), "Elist")
})