context("reduce methods")

test_that("reduce orders columns of each experiment the same", {

  check_reduce_colorder <- function(mae){
    maeredux = reduce(mae)
    maeassay = assay(maeredux)
    matchednames <- lapply(maeassay, function(x){
      sampleMap(maeredux)[match(colnames(x), sampleMap(maeredux)$colname), "primary"]
    })
    for (i in seq_along(2:length(matchednames))){
      testthat::expect_equal(matchednames[[i]], matchednames[[1]])
    }
  }

  example("MultiAssayExperiment")
  experiments(myMultiAssayExperiment)[[1]] <- experiments(myMultiAssayExperiment)[[1]][, 4:1]
  check_reduce_colorder(myMultiAssayExperiment)
  check_reduce_colorder(myMultiAssayExperiment[, , 1:2])
  check_reduce_colorder(myMultiAssayExperiment[, , 1:3])


})
