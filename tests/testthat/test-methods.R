context("Methods tests for methods checks")

library(SummarizedExperiment)
setClass("MAEshell", slots = c(a = "RangedSummarizedExperiment"))
setMethod("assay", "MAEshell", function(x) {assay(x@a)})
setMethod("rownames", "MAEshell", function(x) {rownames(x@a)})
setMethod("colnames", "MAEshell", function(x) {colnames(x@a)})
setMethod("[", c("MAEshell", "ANY", "ANY"), function(x, i, j) {x@a[i, j]})
setMethod("dim", "MAEshell", function(x) {c(length(x@a), length(assay(x@a)[1,]))})

nrows <- 5; ncols <- 4
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rowRanges <- GRanges(rep(c("chr1", "chr2"), c(2, nrows - 2)),
                     IRanges(floor(runif(nrows, 1e5, 1e6)), width=100),
                     strand=sample(c("+", "-"), nrows, TRUE),
                     feature_id=sprintf("ID\\%03d", 1:nrows))
names(rowRanges) <- letters[1:5]
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 2),
                     row.names= c("mysnparray1", "mysnparray2",
                                  "mysnparray3", "mysnparray4"))
rse <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            rowRanges=rowRanges, colData=colData)

masPheno <- data.frame(sex=c("M", "F", "M", "F"),
    age=38:41,
    row.names=c("Jack", "Jill", "Bob", "Barbara"))

aShell <- new("MAEshell", a = rse)
newShell <- list(myShell = aShell)
rangemap <-
   DataFrame(primary = c("Jack", "Jill", "Bob", "Barbara"),
              assay = c("mysnparray1", "mysnparray2", "mysnparray3",
                        "mysnparray4"), assayname = Rle("myShell"))

test_that("the methods check out", {
  expect_true(is(experiments(newShell), "ExperimentList"))
  expect_true(is(MultiAssayExperiment(newShell, masPheno, rangemap), "MultiAssayExperiment"))
})
