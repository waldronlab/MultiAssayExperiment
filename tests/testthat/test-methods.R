context("Methods tests for methods checks")

test_that("the methods check out", {
setClass("MAEshell", slots = c(assay = "matrix"))
setMethod("assay", "MAEshell", function(x) {x@assay})
setMethod("dimnames", "MAEshell", function(x) dimnames(x@assay))
setMethod("[", c("MAEshell", "ANY", "ANY"), function(x, i, j, ...) {
    if (!missing(i) && !missing(j))
        x@assay[i, j]
    else if (missing(i))
        x@assay[, j]
    else if (missing(j))
        x@assay[i, ]
    else x
})
setMethod("dim", "MAEshell", function(x) {
    dim(x@assay)
})

nrows <- 5
ncols <- 4
sampNames <- c("mysnparray1", "mysnparray2", "mysnparray3", "mysnparray4")
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows,
    dimnames = list(NULL, sampNames))
masPheno <- data.frame(sex=c("M", "F", "M", "F"),
    age=38:41,
    row.names=c("Jack", "Jill", "Bob", "Barbara"))

aShell <- new("MAEshell", assay = counts)
newShell <- list(myShell = aShell)
rangemap <-
    DataFrame(assay = factor("myShell"),
              primary = c("Jack", "Jill", "Bob", "Barbara"),
              colname = sampNames)

  expect_true(is(ExperimentList(newShell), "ExperimentList"))
  expect_true(is(MultiAssayExperiment(newShell, masPheno, rangemap),
                 "MultiAssayExperiment"))
})
