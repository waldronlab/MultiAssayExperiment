context("Subset tests")

example("MultiAssayExperiment")

test_that("MultiAssayExperiment length remains the same after subset with list",
{
    rowSubsettor <- rownames(myMultiAssayExperiment)[
        LogicalList(c(TRUE, FALSE), c(FALSE, TRUE))]
    subsettor2 <- rownames(myMultiAssayExperiment)

    expect_true(length(rowSubsettor) !=
        length(myMultiAssayExperiment[rowSubsettor, ]))
    expect_equal(myMultiAssayExperiment[subsettor2, ], myMultiAssayExperiment)
})

test_that("subsetByRow works with lists", {
    rowSubsettor <- rownames(myMultiAssayExperiment)[
        LogicalList(c(TRUE, FALSE), c(FALSE, TRUE))]
    noAffy <- list(noAffy = 1:5)
    expect_error(subsetByRow(experiments(myMultiAssayExperiment), noAffy))
    ## list-like subsets will preserve the original length of the object
    expect_equal(length(subsetByRow(myMultiAssayExperiment, rowSubsettor)),
        length(myMultiAssayExperiment))
})

test_that("assay subsets work", {
    noAffy <- list(noAffy = 1:5)
    expect_error(experiments(myMultiAssayExperiment)[noAffy])
    expect_error(subsetByAssay(myMultiAssayExperiment, noAffy))
    expect_equal(length(subsetByAssay(myMultiAssayExperiment, "Affy")),
        length("Affy"))
    ## check order
    expect_identical(names(subsetByAssay(myMultiAssayExperiment,
        rev(names(myMultiAssayExperiment)))),
        rev(names(myMultiAssayExperiment)))
})

test_that("drop argument works", {
    colList1 <- colnames(myMultiAssayExperiment)
    colList1[[2L]] <- character(0L)

    colList2 <- colList1
    colList2[[4L]] <- character(0L)

    rowList1 <- rownames(myMultiAssayExperiment)
    rowList1[[3L]] <- character(0L)

    rowList2 <- rowList1
    rowList2[[1L]] <- character(0L)

    fullLength <- length(myMultiAssayExperiment)
    minusOne <- length(myMultiAssayExperiment) - 1L
    minusTwo <- length(myMultiAssayExperiment) - 2L
    expect_equal(
        length(myMultiAssayExperiment[, colList1, drop = TRUE]), minusOne
    )
    expect_equal(
        length(myMultiAssayExperiment[, colList2, drop = TRUE]), minusTwo
    )
    expect_equal(
        length(myMultiAssayExperiment[, colList1, drop = FALSE]), fullLength
    )
    expect_equal(
        length(myMultiAssayExperiment[, colList2, drop = FALSE]), fullLength
    )
    expect_equal(
        length(myMultiAssayExperiment[rowList1, drop = TRUE]), minusOne
    )
    expect_equal(
        length(myMultiAssayExperiment[rowList2, drop = TRUE]), minusTwo
    )
    expect_equal(
        length(myMultiAssayExperiment[rowList1, drop = FALSE]), fullLength
    )
    expect_equal(length(myMultiAssayExperiment[FALSE, drop = TRUE]), 0L)
    expect_equal(
        length(myMultiAssayExperiment[FALSE, drop = FALSE]), fullLength
    )
    expect_equal(length(myMultiAssayExperiment[, FALSE, drop = TRUE]), 0L)
    expect_equal(
        length(myMultiAssayExperiment[, FALSE, drop = FALSE]), fullLength
    )
})

test_that("subsetByColumn works with lists", {
    affySub <- list(Affy = 1:2)
    affySimple <- List(affySub)
    expect_equal(length(myMultiAssayExperiment[, affySub, ]), length(affySub))
    expect_equal(length(myMultiAssayExperiment[, affySimple, ]),
                 length(affySimple))
    expect_equal(length(experiments(myMultiAssayExperiment)[affySub]),
                 length(affySub))
    expect_equal(length(experiments(myMultiAssayExperiment)[affySimple]),
                 length(affySimple))
})
test_that("subsetByColumn works with incomplete lists", {
    fuLL <- colnames(myMultiAssayExperiment)
    cLIST <- fuLL[c(1,3)]

    expect_equal(length(myMultiAssayExperiment[, cLIST, , drop = FALSE]),
        length(fuLL))
    expect_equal(length(myMultiAssayExperiment[, cLIST, , drop = TRUE]),
        length(cLIST))
})

test_that("subsetByColData works as intended", {
    trues <- sum(myMultiAssayExperiment$sex == "M")
    expect_equal(
        nrow(colData(subsetByColData(myMultiAssayExperiment,
            myMultiAssayExperiment$sex == "M"))),
        trues
    )
    expect_equal(
        nrow(colData(myMultiAssayExperiment[, myMultiAssayExperiment$sex == "M"])),
        trues
    )
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
    expectedOrder <- IRanges::CharacterList(
        A = c("row1", "row2"),
        B = c("row1", "row2"))
    expect_identical(rownames(mae[c("row1", "row2")]), expectedOrder)
})
