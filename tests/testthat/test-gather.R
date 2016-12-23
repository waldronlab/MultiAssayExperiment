context("gather methods")

example("MultiAssayExperiment")

test_that("gather returns a DataFrame", {

    testMat <- matrix(seq_len(20), nrow = 4, ncol = 5, byrow = TRUE,
                      dimnames = list(LETTERS[1:4], letters[1:5]))

    testESet <- Biobase::ExpressionSet(testMat)

    matDF <- gather(testMat)
    ESetDF <- gather(testESet)
    RRADF <- gather(myRRA)

    expect_true(is(matDF, "DataFrame"))
    expect_true(is(ESetDF, "DataFrame"))
    expect_true(is(RRADF, "DataFrame"))
})
