context("helper methods")

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
