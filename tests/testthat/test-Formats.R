context("long and wide format methods")

test_that("longFormat returns a DataFrame", {
    testMat <- matrix(seq_len(20), nrow = 4, ncol = 5, byrow = TRUE,
                      dimnames = list(LETTERS[1:4], letters[1:5]))

    testESet <- Biobase::ExpressionSet(testMat)

    matDF <- longFormat(testMat)
    ESetDF <- longFormat(testESet)

    expect_true(is(matDF, "data.frame"))
    expect_true(is(ESetDF, "data.frame"))
})

test_that("longFormat returns specified colData column and proper dimensions", {
    example("MultiAssayExperiment")

    longDF <- longFormat(myMultiAssayExperiment, colDataCols = "sex")
    wideDF <- wideFormat(myMultiAssayExperiment, colDataCols = "sex")

    expect_true("sex" %in% names(longDF))
    expect_true("sex" %in% names(wideDF))
    expect_equal(nrow(wideDF), nrow(colData(myMultiAssayExperiment)))
})

test_that("wideFormat returns primary column order identical to colData rownames", {
    data("miniACC")
    acc <- miniACC["EZH2", , ]
    wideacc <- wideFormat(acc)
    expect_true(identical(rownames(colData(acc)),
        as.character(wideacc[["primary"]])))
})
