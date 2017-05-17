context("Subset tests")

example("MultiAssayExperiment")

test_that("subsettor length is of the same as MultiAssayExperiment", {
    rowSubsettor <- lapply(rownames(myMultiAssayExperiment)[1:2],
                           function(a) { sample(a, 1) })
    subsettor2 <- rownames(myMultiAssayExperiment)

    expect_false(length(rowSubsettor) == length(myMultiAssayExperiment))
    expect_error(myMultiAssayExperiment[rowSubsettor, ])
    expect_equal(myMultiAssayExperiment[subsettor2, ], myMultiAssayExperiment)
})

test_that("drop argument works", {
    colList <- colnames(myMultiAssayExperiment)
    colList[[2L]] <- character(0L)

    rowList <- rownames(myMultiAssayExperiment)
    rowList[[3L]] <- character(0L)

    fullLength <- length(myMultiAssayExperiment)
    minusOne <- length(myMultiAssayExperiment)-1L
    expect_equal(length(myMultiAssayExperiment[, colList, drop = TRUE]),
                 minusOne)
    expect_equal(length(myMultiAssayExperiment[, colList, drop = FALSE]),
                 fullLength)
    expect_equal(length(myMultiAssayExperiment[rowList, drop = TRUE]),
                 minusOne)
    expect_equal(length(myMultiAssayExperiment[rowList, drop = FALSE]),
                 fullLength)
    expect_equal(length(myMultiAssayExperiment[FALSE, drop = TRUE]), 0L)
    expect_equal(length(myMultiAssayExperiment[FALSE, drop = FALSE]),
                 fullLength)
    expect_equal(length(myMultiAssayExperiment[, FALSE, drop = TRUE]), 0L)
    expect_equal(length(myMultiAssayExperiment[, FALSE, drop = FALSE]),
                 fullLength)
})

test_that("subsetByColumns works with lists", {
    affySub <- list(Affy = 1:2)
    affySimple <- List(affySub)
    expect_equal(length(myMultiAssayExperiment[, affySub, ]), length(affySub))
    expect_equal(length(myMultiAssayExperiment[, affySimple, ]),
                 length(affySimple))
})

test_that("subsetByColData works as intended", {
    trues <- sum(myMultiAssayExperiment$sex == "M")
    expect_equal(nrow(colData(subsetByColData(myMultiAssayExperiment,
                                          myMultiAssayExperiment$sex == "M"))),
                 trues)
    expect_equal(nrow(colData(
        myMultiAssayExperiment[, myMultiAssayExperiment$sex == "M"])), trues)
})

test_that("MultiAssayExperiment subsetting works with NULL rownames", {
    nrows <- 200
    ncols <- 6
    counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
    rowRanges <- GRanges(rep(c("chr1", "chr2"), c(50, 150)),
                         IRanges(floor(runif(200, 1e5, 1e6)), width=100),
                         strand=sample(c("+", "-"), 200, TRUE),
                         feature_id=sprintf("ID%03d", 1:200))
    colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                         row.names=LETTERS[1:6])
    rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                                rowRanges=rowRanges, colData=colData)
    maeRSE <- MultiAssayExperiment(list(rag = rse))
    expect_true(is.character(rownames(maeRSE)[[1L]]))
    expect_true(length(rownames(maeRSE)[[1L]]) == 0L)
})

test_that("subsetting by row keeps the order of the subsettor", {
    mdat1 <- matrix(c(1, 11, 2, 12, 3, 13), nrow = 2,
        dimnames = list(c("row1", "row2"), c("C1", "C2", "C3")))
    mdat2 <- mdat1[2:1, ]
    mae <- MultiAssayExperiment(list(A = mdat1, B = mdat2))
    expectedOrder <- IRanges::CharacterList(A = c("row1", "row2"),
                                            B = c("row1", "row2"))
    expect_identical(rownames(mae[c("row1", "row2")]), expectedOrder)
})
