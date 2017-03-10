test_that("combine c function works", {
    mat <- matrix(rnorm(20), ncol = 4, dimnames = list(NULL, LETTERS[1:4]))
    myMulti <- MultiAssayExperiment(list(express = mat))
    mat2 <- matrix(rnorm(20), ncol = 4, dimnames = list(NULL, letters[1:4]))
    expect_error(c(myMulti, mat2)) # Unnamed list element
    expect_error(c(myMulti, express2 = mat2)) # no sampleMap or mapFrom arg
    expect_true(validObject(c(myMulti, express2 = mat2, mapFrom = 1L)))
    expect_true(is(c(myMulti, express2 = mat2, mapFrom = 1L),
                   "MultiAssayExperiment"))
})
