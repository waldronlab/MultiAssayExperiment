test_that("MultiAssayExperiment constructors work", {
    expect_true(validObject(MultiAssayExperiment()))

    expect_error(MultiAssayExperiment(list(matrix(0, 0, 0))),
                 label="unnamed ExperimentList")
})
