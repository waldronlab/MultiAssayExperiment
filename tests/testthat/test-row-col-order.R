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
            expect_equal(matchednames[[i]], matchednames[[1]])
    }

    mae[[1]] <- mae[[1]][, 4:1]
    check_mergeReps_colorder(mae)
    check_mergeReps_colorder(mae[, , 1:2])
    check_mergeReps_colorder(mae[, , 1:3])

})

test_that("row order identical after matching", {
    introws <- Reduce(intersect, rownames(mae))
    intMultiRow <- rownames(intersectRows(mae))
    expect_true(
        all(vapply(intMultiRow,
            function(assayDat) {
                identical(assayDat, introws)
            }, logical(1L))
        )
    )

    reord <- mendoapply(
        FUN = rev,
        rownames(mae)[rep(list(c(TRUE, FALSE)), 4)]
    )
    mord <- mae[reord, , drop = FALSE]
    expect_true(
        all(mapply(identical, rownames(mord), reord))
    )
})
