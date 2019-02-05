context("Endomorphic tests of builder functions")

test_that("builders return empty", {
    expect_true(isEmpty(MultiAssayExperiment()))
    expect_true(isEmpty(ExperimentList()))
})

test_that("builders return appropriate class", {
    expect_true(is(MultiAssayExperiment(), "MultiAssayExperiment"))
})

test_that("slots are of appropriate class", {
    expect_true(is(colData(MultiAssayExperiment()), "DataFrame"))
    expect_true(is(sampleMap(MultiAssayExperiment()), "DataFrame"))
    expect_true(is(ExperimentList(), "ExperimentList"))
    expect_true(is(metadata(MultiAssayExperiment()), "list"))
})
