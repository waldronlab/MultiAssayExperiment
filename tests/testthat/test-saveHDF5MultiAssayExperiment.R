context("check saveHDF5MultiAssayExperiment")

test_that("saveHDF5MultiAssayExperiment is working", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]

    testDir <- file.path(tempdir(), "test_mae")
    dir.create(testDir)
    if (!file.exists(file.path(testDir, "miniACC_experiments.h5")))
        saveHDF5MultiAssayExperiment(miniACC, dir = testDir, replace = TRUE)

    expect_identical(
        list.files(testDir), c("miniACC_experiments.h5", "miniACC_mae.rds")
    )

    x <- loadHDF5MultiAssayExperiment(dir = testDir)
    expect_true(is(x, "MultiAssayExperiment"))
    vapply(experiments(x), function(expr) {
        expect_true(is(expr, "HDF5Matrix"))
    }, logical(1L))

    testDir <- file.path(tempdir(), "test_mae")
    saveHDF5MultiAssayExperiment(
        miniACC, dir = testDir, prefix = "", replace = TRUE
    )
    expect_identical(
        list.files(testDir), c("experiments.h5", "mae.rds")
    )
})

test_that("array assays work with saveHDF5MultiAssayExperiment", {
    A <- array(1:24, 4:2, dimnames =
        list(letters[1:4], LETTERS[1:3], c("A", "B"))
    )
    B <- matrix(1:12, ncol=3, dimnames = list(LETTERS[1:4], letters[1:3]))
    se2 <- SummarizedExperiment(list(A=A, B=B))
    mae2 <- MultiAssayExperiment(ExperimentList(one_more_se=se2))
    saveHDF5MultiAssayExperiment(mae2, "mae2", replace = TRUE)
    loadHDF5MultiAssayExperiment("mae2")
})
