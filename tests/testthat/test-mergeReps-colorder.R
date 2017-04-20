context("helper methods")

test_that("merging replicates orders columns of each experiment the same", {

check_mergeReps_colorder <- function(mae) {
    maeredux = mergeReplicates(intersectColumns(mae))
    maeassay = assays(maeredux)
    matchednames <- lapply(maeassay, function(x) {
      sampleMap(maeredux)[match(colnames(x),
                                sampleMap(maeredux)$colname),
                          "primary"]
    })
    for (i in seq_along(2:length(matchednames)))
      testthat::expect_equal(matchednames[[i]], matchednames[[1]])
}

    example("MultiAssayExperiment")
    myMultiAssayExperiment[[1]] <- myMultiAssayExperiment[[1]][, 4:1]
    check_mergeReps_colorder(myMultiAssayExperiment)
    check_mergeReps_colorder(myMultiAssayExperiment[, , 1:2])
    check_mergeReps_colorder(myMultiAssayExperiment[, , 1:3])

})
