context("Subset tests")

example("MultiAssayExperiment")

colList <- colnames(myMultiAssayExperiment)
colList[[2L]] <- character(0L)

rowList <- rownames(myMultiAssayExperiment)
rowList[[3L]] <- character(0L)

rowSubsettor <- lapply(rownames(myMultiAssayExperiment)[1:2], function(a) {
  sample(a, 1)
})
subsettor2 <- rownames(myMultiAssayExperiment)

test_that("subsettor length is of the same as MultiAssayExperiment", {
  expect_false(length(rowSubsettor) == length(myMultiAssayExperiment))
  expect_error(myMultiAssayExperiment[rowSubsettor, ])
  expect_equal(myMultiAssayExperiment[subsettor2, ], myMultiAssayExperiment)
})

test_that("drop argument works", {
  expect_equal(length(myMultiAssayExperiment[, colList, drop = TRUE]), 3L)
  expect_equal(length(myMultiAssayExperiment[, colList, drop = FALSE]), 4L)
  expect_equal(length(myMultiAssayExperiment[rowList, drop = TRUE]), 3L)
  expect_equal(length(myMultiAssayExperiment[rowList, drop = FALSE]), 4L)
  expect_equal(length(myMultiAssayExperiment[FALSE, drop = TRUE]), 0L)
  expect_equal(length(myMultiAssayExperiment[FALSE, drop = FALSE]), 4L)
  expect_equal(length(myMultiAssayExperiment[, FALSE, drop = TRUE]), 0L)
  expect_equal(length(myMultiAssayExperiment[, FALSE, drop = FALSE]), 4L)
})

test_that("subsetByColumns works with lists", {
    affySub <- list(Affy = 1:2)
    affySimple <- List(affySub)
    expect_equal(length(myMultiAssayExperiment[, affySub, ]), length(affySub))
    expect_equal(length(myMultiAssayExperiment[, affySimple, ]),
                 length(affySimple))
})

test_that("subsetBypData works as intended", {
    trues <- sum(myMultiAssayExperiment$sex == "M")
    expect_equal(nrow(pData(subsetBypData(myMultiAssayExperiment,
                                          myMultiAssayExperiment$sex == "M"))),
                 trues)
    expect_equal(nrow(pData(
        myMultiAssayExperiment[, myMultiAssayExperiment$sex == "M"])), trues)
})
