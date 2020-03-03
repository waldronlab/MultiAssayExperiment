context("Checking for consistency in colnames and rownames")

test_that("rownames and colnames return CompressedCharacterList", {
    expect_true(is(rownames(mae), "CompressedCharacterList"))
    expect_true(is(colnames(mae), "CompressedCharacterList"))
})
