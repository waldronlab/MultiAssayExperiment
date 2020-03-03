context("drop argument for rows works as it should")

test_that("when drop is TRUE, length of ExperimentList is correct", {
    expect_true(length(mae[FALSE, , , drop = TRUE]) == 0L)
})
test_that("when drop is FALSE, length is the same", {
    expect_true(length(mae[FALSE, , , drop = FALSE]) == length(mae))
})

context("drop arugment for columns works as it should")

test_that("when drop is TRUE, length of ExperimentList is correct", {
    expect_true(length(mae[, FALSE, , drop = TRUE]) == 0L)
})
test_that("when drop is FALSE, length is the same", {
    expect_true(length(mae[, FALSE, , , drop = FALSE]) == length(mae))
})
