context("rearrange methods")

example("MultiAssayExperiment")

test_that("rearrange returns a DataFrame", {

    testMat <- matrix(seq_len(20), nrow = 4, ncol = 5, byrow = TRUE,
                      dimnames = list(LETTERS[1:4], letters[1:5]))

    testESet <- Biobase::ExpressionSet(testMat)

    matDF <- rearrange(testMat)
    ESetDF <- rearrange(testESet)
    RRADF <- rearrange(myRRA)

    expect_true(is(matDF, "DataFrame"))
    expect_true(is(ESetDF, "DataFrame"))
    expect_true(is(RRADF, "DataFrame"))
})
