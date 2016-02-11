context("Endomorphic tests of builder functions")

test_that("builders return empty", {
  expect_true(isEmpty(MultiAssayExperiment()))
  expect_true(isEmpty(RangedRaggedAssay()))
  expect_true(isEmpty(Elist()))
})

test_that("builders return appropriate class", {
  expect_true(is(MultiAssayExperiment(), "MultiAssayExperiment"))
  expect_true(is(RangedRaggedAssay(), "RangedRaggedAssay"))
})

test_that("slots are of appropriate class", {
  expect_true(is(pData(MultiAssayExperiment()), "DataFrame"))
  expect_true(is(sampleMap(MultiAssayExperiment()), "DataFrame"))
  expect_true(is(Elist(), "Elist"))
  expect_null(metadata(MultiAssayExperiment()))
})
