context("c function expands ExperimentList and modifies sampleMap")

example("MultiAssayExperiment")

## BiocInstaller::biocLite("Bioconductor/RaggedExperiment")

library(SummarizedExperiment)
library(RaggedExperiment)
example("SummarizedExperiment")
example("RaggedExperiment")

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

test_that("combine c function works on multiple objects", {
    rse1 <- rse1[, 1:4]
    addMap <- DataFrame(assay = "RExp",
                        primary = c("Jack", "Jill"),
                        colname = c("sample1", "sample2"))
    addMap <- rbind(addMap, DataFrame(assay = "myRSE",
                                      primary = c("Jack", "Jill",
                                                  "Bob", "Barbara"),
                                      colname = LETTERS[1:4]))
    expect_true(validObject(c(myMultiAssayExperiment, RExp = re3,
                              myRSE = rse1, sampleMap = addMap)))
})

