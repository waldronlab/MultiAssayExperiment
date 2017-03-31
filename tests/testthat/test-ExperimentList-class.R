context("ExperimentList tests")

test_that("ExperimentList constructors work", {
    expect_identical(names(ExperimentList()), character())
    expect_true(validObject(ExperimentList()))

    expect_true(validObject(ExperimentList(list(m=matrix(0, 0, 0)))))
    expect_error(ExperimentList(list(matrix(0, 0, 0))),
                 label="unnamed list element")

    m <- matrix(0, 2, 2, dimnames=list(letters[1:2], letters[1:2]))
    expect_true(validObject(ExperimentList(list(m=m))))
})
