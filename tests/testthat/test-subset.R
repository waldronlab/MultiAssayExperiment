context("Subset tests")

test_that("MultiAssayExperiment length remains the same after subset with list",
{
    rowindex <- rownames(mae)[LogicalList(c(TRUE, FALSE), c(FALSE, TRUE))]
    index2 <- rownames(mae)

    ## drop = TRUE by default
    expect_identical(
        length(mae), length(mae[rowindex, drop = TRUE])
    )
    expect_equal(rownames(mae[index2, ]), index2)
})

test_that("subsetByRow works with lists", {
    fullROWS <- rownames(mae)
    rows <- fullROWS[LogicalList(c(TRUE, FALSE), c(FALSE, TRUE))]
    noAffy <- list(noAffy = 1:5)
    expect_error(subsetByRow(mae, noAffy))
    expect_error(subsetByRow(experiments(mae), noAffy))
    ## list-like subsets will preserve the original length of the object
    expect_equal(length(subsetByRow(mae, rows)), length(mae))
    expect_equal(length(mae[rows, ]), length(mae))
})

test_that("subsetByRow keeps assay order in MultiAssayExperiment", {
    rows <- rownames(mae)
    ## re-arrange experiment rownames
    rows <- rows[c(2, 3, 1, 4)]
    newrows <- subsetByRow(mae, rows)
    expect_identical(names(mae), names(newrows))
    expect_identical(names(mae), names(mae[rows, ]))
})

test_that("assay subsets work", {
    noAffy <- list(noAffy = 1:5)
    expect_error(experiments(mae)[noAffy])
    expect_error(subsetByAssay(mae, noAffy))
    expect_equal(length(subsetByAssay(mae, "Affy")),
        length("Affy"))
    ## check order
    expect_identical(names(subsetByAssay(mae, rev(names(mae)))),
        rev(names(mae)))
})

test_that("drop argument works", {
    colList1 <- colnames(mae)
    colList1[[2L]] <- character(0L)

    colList2 <- colList1
    colList2[[4L]] <- character(0L)

    rowList1 <- rownames(mae)
    rowList1[[3L]] <- character(0L)

    rowList2 <- rowList1
    rowList2[[1L]] <- character(0L)

    fullLength <- length(mae)
    minusOne <- length(mae) - 1L
    minusTwo <- length(mae) - 2L
    expect_equal(
        length(mae[, colList1, drop = TRUE]), minusOne
    )
    expect_equal(
        length(mae[, colList2, drop = TRUE]), minusTwo
    )
    expect_equal(
        length(mae[, colList1, drop = FALSE]), fullLength
    )
    expect_equal(
        length(mae[, colList2, drop = FALSE]), fullLength
    )
    expect_equal(
        length(mae[rowList1, drop = TRUE]), minusOne
    )
    expect_equal(
        length(mae[rowList2, drop = TRUE]), minusTwo
    )
    expect_equal(
        length(mae[rowList1, drop = FALSE]), fullLength
    )
    expect_equal(length(mae[FALSE, drop = TRUE]), 0L)
    expect_equal(
        length(mae[FALSE, drop = FALSE]), fullLength
    )
    expect_equal(length(mae[, FALSE, drop = TRUE]), 0L)
    expect_equal(
        length(mae[, FALSE, drop = FALSE]), fullLength
    )
})

test_that("subsetByColumn keeps order in index", {
    cols <- colnames(mae)
    cols <- cols[c(2, 3, 1, 4)]
    newcols <- subsetByColumn(mae, cols)
    expect_identical(names(mae), names(newcols))
    expect_identical(names(mae), names(mae[, cols]))
})

test_that("subsetByColumn works with lists", {
    affySub <- list(Affy = 1:2)
    affySimple <- List(affySub)

    expect_equal(length(mae[, affySub, ]), length(mae))
    expect_equal(length(mae[, affySimple, ]), length(mae))
    expect_equal(length(experiments(mae)[affySub]), length(affySub))
    expect_equal(length(experiments(mae)[affySimple]), length(affySimple))

    # incomplete lists
    fuLL <- colnames(mae)
    cLIST <- fuLL[c(1,3)]
    pLIST <-
        list(Affy = 1:4, Methyl450k = integer(), RNASeqGene = 1:4, GISTIC = 1:3)

    expect_equal(length(mae[, cLIST, , drop = FALSE]), length(fuLL))
    expect_equal(length(mae[, cLIST, , drop = TRUE]), length(mae))

    ## test zero length subset in list
    expect_equal(length(mae[, pLIST, drop = FALSE]), length(fuLL))
    expect_equal(length(mae[, pLIST, drop = TRUE]), length(fuLL)-1L)
})

test_that("subsetByColData works as intended", {
    mae0 <- mae
    mae0 <- subsetByColData(mae, mae$sex == "M")
    expect_identical(
        colData(mae0),
        colData(mae)[mae$sex == "M", ]
    )

    ## probably a harmonize check..
    newsm <- sampleMap(mae)
    newsm <-
        newsm[newsm$primary %in% rownames(colData(mae)[mae$sex == "M", ]), ]
    expect_identical(
        sampleMap(mae0),
        newsm
    )

    longLogi <- rep(c(TRUE, FALSE), 3)
    expect_error(
        subsetByColData(mae, longLogi)
    )

    chars <- rownames(colData(mae))[c(TRUE, FALSE)]
    mae0 <- subsetByColData(mae, chars)
    expect_identical(
        colData(mae0),
        colData(mae)[chars, ]
    )

    mae0 <- subsetByColData(mae, rev(chars))
    expect_identical(
        colData(mae0),
        colData(mae)[rev(chars), ]
    )

    mae0 <- subsetByColData(mae, c(1, 4))
    expect_identical(
        colData(mae0),
        colData(mae)[c(1, 4), ]
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

test_that("subsetting by row keeps the order of the index", {
    mdat1 <- matrix(c(1, 11, 2, 12, 3, 13), nrow = 2,
        dimnames = list(c("row1", "row2"), c("C1", "C2", "C3")))
    mdat2 <- mdat1[2:1, ]
    mae <- MultiAssayExperiment(list(A = mdat1, B = mdat2))
    expectedOrder <- IRanges::CharacterList(
        A = c("row1", "row2"),
        B = c("row1", "row2"))
    expect_identical(rownames(mae[c("row1", "row2")]), expectedOrder)
})

test_that("replacement methods work", {
    ## test double bracket replacement
    exse <- ExpList[["Affy"]]
    assay(exse) <- log(assay(ExpList[["Affy"]]))
    expect_error(
        mae[[1L]] <- SimpleList(Affy = exse)
    )
    mae[[1L]] <- exse
    expect_identical(assay(mae[[1L]]), assay(exse))

    ## test multiple replacement as list or List
    replace2 <- log(mae[[2L]] + 1)
    mae[, , 1:2] <- SimpleList(exse, replace2)
    expect_identical(mae[[2]], replace2)

    ## test multiple replacements with names
    replaceExList <- ExperimentList(logAffy = exse, logMeth = replace2)
    mae[, , 1:2] <- replaceExList
    expect_identical(mae[[1]], exse)
    expect_identical(mae[[2]], replace2)
    expect_identical(names(mae[, , 1:2]), names(replaceExList))
})
