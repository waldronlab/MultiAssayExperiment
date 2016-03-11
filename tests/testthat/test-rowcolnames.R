context("Checking for consistency in colnames and rownames")

example("MultiAssayExperiment")

test_that("rownames and colnames return CompressedCharacterList", {
    expect_true(is(rownames(myMultiAssayExperiment),
                   "CompressedCharacterList"))
    expect_true(is(colnames(myMultiAssayExperiment),
                   "CompressedCharacterList"))
})
