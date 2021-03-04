context("Class checks for constructor functions")

library(GenomicRanges)
gr1 <-
  GRanges(seqnames = "chr3", ranges = IRanges(58000000, 59502360),
          strand = "+", score = 5L, GC = 0.45)
gr2 <-
  GRanges(seqnames = c("chr3", "chr3"),
          ranges = IRanges(c(58493000, 3), width=9000),
          strand = c("+", "-"), score = 3:4, GC = c(0.3, 0.5))
gr3 <-
  GRanges(seqnames = c("chr1", "chr2"),
          ranges = IRanges(c(1, 4), c(3, 9)),
          strand = c("-", "-"), score = c(6L, 2L), GC = c(0.4, 0.1))

grl <- GRangesList("gr1" = gr1, "gr2" = gr2, "gr3" = gr3)
names(grl) <- c("snparray1", "snparray2", "snparray3")
RagExp <- RaggedExperiment::RaggedExperiment(grl)

arraydat <- matrix(seq(101, 108), ncol=4,
                   dimnames = list(
                     c("ENST00000294241", "ENST00000355076"),
                     c("array1", "array2", "array3", "array4")
                   ))
arraypdat <- as(data.frame(
  slope53 = rnorm(4),
  row.names = c("array1", "array2", "array3", "array4")),
  "AnnotatedDataFrame")
exprdat <- Biobase::ExpressionSet(assayData=arraydat, phenoData=arraypdat)

assayList <- list(RagExp, exprdat)
names(assayList) <- c("CNVgistic", "Affy")
ExpList <- ExperimentList(assayList)

test_that("the appropriate class is returned", {
  expect_true(is(ExpList, "ExperimentList"))
})
