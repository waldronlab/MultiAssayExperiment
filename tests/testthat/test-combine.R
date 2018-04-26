context("c function expands ExperimentList and modifies sampleMap")

test_that("combine c function works", {
    mat <- matrix(rnorm(20), ncol = 4, dimnames = list(NULL, LETTERS[1:4]))
    myMulti <- MultiAssayExperiment(list(express = mat))
    mat2 <- matrix(rnorm(20), ncol = 4, dimnames = list(NULL, letters[1:4]))
    expect_error(c(myMulti, mat2)) # Unnamed list element
    # no sampleMap or mapFrom arg
    expect_error(c(myMulti, express2 = mat2))
    expect_true(validObject(c(myMulti, express2 = mat2, mapFrom = 1L)))
    expect_true(is(c(myMulti, express2 = mat2, mapFrom = 1L),
        "MultiAssayExperiment"))
})

test_that("combine c function works on multiple objects", {
    example("MultiAssayExperiment")

    library(SummarizedExperiment)
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
    rse <- rse[, 1:4]

    library(RaggedExperiment)
    sample1 <- GRanges(c(a = "chr1:1-10:-", b = "chr1:11-18:+"), score = 1:2)
    sample2 <- GRanges(c(c = "chr2:1-10:-", d = "chr2:11-18:+"), score = 3:4)
    colDat <- DataFrame(id = 1:2)
    grl <- GRangesList(sample1 = sample1, sample2 = sample2)
    re <- RaggedExperiment(grl, colData = colDat)

    addMap <- DataFrame(assay = "RExp",
        primary = c("Jack", "Jill"),
        colname = c("sample1", "sample2"))
    addMap <- rbind(addMap, DataFrame(assay = "myRSE",
        primary = c("Jack", "Jill", "Bob", "Barbara"),
        colname = LETTERS[1:4]))
    expect_true(validObject(c(myMultiAssayExperiment, RExp = re,
        myRSE = rse, sampleMap = addMap)))
})

