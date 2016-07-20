test_that("MultiAssayExperiment constructors work", {
    expect_true(validObject(MultiAssayExperiment()))

    expect_error(MultiAssayExperiment(list(matrix(0, 0, 0))),
                 label="unnamed ExperimentList")

    obs <- MultiAssayExperiment(list(m=matrix(0, 0, 0)))
    expect_true(validObject(obs))

    obs <- MultiAssayExperiment(list(stack=GenomicFiles::VcfStack()))
    expect_true(validObject(obs))
})
