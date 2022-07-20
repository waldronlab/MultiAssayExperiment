context("check saveHDF5MultiAssayExperiment")

test_that("saveHDF5MultiAssayExperiment is working", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]

    testDir <- file.path(tempdir(), "test_mae")
    saveHDF5MultiAssayExperiment(miniACC, dir = testDir, replace = TRUE)

    expect_identical(
        list.files(testDir), c("miniACC_experiments.h5", "miniACC_mae.rds")
    )

    x <- loadHDF5MultiAssayExperiment(dir = testDir)
    expect_true(validObject(x))
    expect_true(is(x, "MultiAssayExperiment"))

    ## test that classes are largely the same (except matrix > HDF5Matrix)
    resclasses <- c(
        rep("SummarizedExperiment", 3), "HDF5Matrix", "SummarizedExperiment"
    )
    names(resclasses) <- names(x)
    expect_equal(
        vapply(experiments(x), class, character(1L)),
        resclasses
    )
    on.exit(unlink(testDir, recursive = TRUE))
})

test_that("prefix argument works as intended", {
    testDir <- file.path(tempdir(), "test_mae")
    saveHDF5MultiAssayExperiment(
        miniACC, dir = testDir, prefix = "", replace = TRUE
    )
    exp_files <- c("experiments.h5", "mae.rds")

    expect_identical(list.files(testDir), exp_files)

    mae <- miniACC
    testDir0 <- file.path(tempdir(), "test_mae0", .Platform$file.sep)
    saveHDF5MultiAssayExperiment(
        mae, dir = testDir0, replace = TRUE
    )
    file.copy(from = file.path(testDir, exp_files), to = testDir0)
    expect_error(loadHDF5MultiAssayExperiment(testDir0))
    expect_true(
        validObject(loadHDF5MultiAssayExperiment(testDir0, "mae"))
    )
    expect_true(
        validObject(loadHDF5MultiAssayExperiment(testDir0, ""))
    )
    on.exit({
        unlink(testDir0, recursive = TRUE)
        unlink(testDir, recursive = TRUE)
    })
})

test_that("array assays work with saveHDF5MultiAssayExperiment", {
    A <- array(1:24, 4:2, dimnames =
        list(letters[1:4], LETTERS[1:3], c("A", "B"))
    )
    B <- matrix(1:12, ncol=3, dimnames = list(letters[1:4], LETTERS[1:3]))
    se2 <- SummarizedExperiment(list(A=A, B=B))
    mae2 <- MultiAssayExperiment(ExperimentList(one_more_se=se2))
    testDir <- file.path(tempdir(), "test_mae")
    expect_warning(
        saveHDF5MultiAssayExperiment(mae2, testDir, replace = TRUE)
    )
    expect_true(
        validObject(loadHDF5MultiAssayExperiment(testDir))
    )
    expect_true(
        validObject(loadHDF5MultiAssayExperiment(testDir, "mae2"))
    )
    on.exit(unlink(testDir, recursive = TRUE))
})

context("check loadHDF5MultiAssayExperiment")

test_that("loadHDF5MultiAssayExperiment is working", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]

    testDir <- file.path(tempdir(), "test_mae")
    on.exit(unlink(testDir, recursive = TRUE))
    
    saveHDF5MultiAssayExperiment(
        miniACC, prefix = "", dir = testDir, replace = TRUE
    )

    expect_true(
        validObject(loadHDF5MultiAssayExperiment(dir = testDir, prefix = ""))
    )
    expect_true(
        validObject(loadHDF5MultiAssayExperiment(dir = testDir))
    )
})

test_that("loadHDF5MultiAssayExperiment prefix input is consistent", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]

    testDir <- file.path(tempdir(), "test_mae")
    on.exit(unlink(testDir, recursive = TRUE))
    
    saveHDF5MultiAssayExperiment(
        miniACC, prefix = "test", dir = testDir, replace = TRUE
    )

    expect_true(
        validObject(
            loadHDF5MultiAssayExperiment(dir = testDir, prefix = "test")
        )
    )
    expect_true(
        validObject(
            loadHDF5MultiAssayExperiment(dir = testDir, prefix = "test_")
        )
    )
    expect_error(
        validObject(
            loadHDF5MultiAssayExperiment(dir = testDir, prefix = "error")
        )
    )
})