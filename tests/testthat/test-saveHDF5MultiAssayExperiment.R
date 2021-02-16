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
