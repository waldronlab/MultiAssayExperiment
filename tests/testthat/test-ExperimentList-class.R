test_that("ExperimentList constructors work", {
    expect_identical(names(ExperimentList()), character())
    expect_true(validObject(ExperimentList()))
})
