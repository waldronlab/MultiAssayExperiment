context("ExperimentList tests")

test_that("ExperimentList constructors work", {
    expect_identical(names(ExperimentList()), character())
    expect_true(validObject(ExperimentList()))

    expect_true(validObject(ExperimentList(list(m=matrix(0, 0, 0)))))
    expect_error(ExperimentList(list(matrix(0, 0, 0))),
                 label="unnamed list element")

    m <- matrix(0, 2, 2, dimnames=list(letters[1:2], letters[1:2]))
    expect_true(validObject(ExperimentList(list(m=m))))
})

test_that("ExperimentList constructor preserves metadata", {
    mcoldf <- DataFrame(AssayNumber=seq_len(length(ExpList)),
        row.names = names(ExpList))
    mcols(ExpList) <- mcoldf

    metalist <- list(Shiny = "Blue Jeans", Old = "Metadata")
    metadata(ExpList) <- metalist

    nexp <- ExperimentList(ExpList)

    expect_identical(mcols(nexp), mcoldf)
    expect_identical(metadata(nexp), metalist)
})

test_that("Metadata is kept in ExperimentList when replacing", {
    ## add metadata columns and metadata
    mcoldf <- DataFrame(AssayNumber=seq_len(length(ExpList)),
        row.names = names(ExpList))
    mcols(ExpList) <- mcoldf

    metalist <- list(Shiny = "Blue Jeans", Old = "Metadata")
    metadata(ExpList) <- metalist

    mae0 <- mae
    experiments(mae0) <- ExpList

    expect_identical(mcols(experiments(mae0)), mcoldf)
    expect_identical(metadata(experiments(mae0)), metalist)
})

test_that("ExperimentList subsetting works", {

    ## subset by rows
    rsub <- rownames(ExpList)[list(1, 2)]
    res <- ExpList[rsub, , ]
    expect_identical(
        rownames(res)[1:2],
        rsub
    )
    rlsub <- S4Vectors::endoapply(rownames(ExpList), `[`, c(TRUE, FALSE))
    res <- ExpList[rlsub, , ]
    expect_identical(
        rownames(res),
        rlsub
    )

    res <- ExpList[1:2, , ]
    expect_identical(
        vapply(res, nrow, integer(1L)),
        stats::setNames(rep(2L, length(res)), names(res))
    )

    res <- ExpList[c(TRUE, FALSE), , ]
    expect_identical(
        vapply(res, nrow, integer(1L)),
        vapply(
            ExpList,
            function(x) {
                sum(!seq_len(nrow(x)) %% 2 == 0L)
            },
            integer(1L)
        )
    )

    ## subset by cols
    csub <- colnames(ExpList)[list(1, 2)]
    res <- ExpList[, csub]
    expect_identical(
        colnames(res)[1:2],
        csub
    )

    clsub <- S4Vectors::endoapply(colnames(ExpList), `[`, c(TRUE, FALSE))
    res <- ExpList[, clsub, ]
    expect_identical(
        colnames(res),
        clsub
    )
    clsub <- S4Vectors::endoapply(colnames(ExpList), `[`, c(TRUE, FALSE))

    ## subset by assay
    asub <- names(ExpList)[c(1,4)]
    res <- ExpList[, , asub]
    expect_identical(
        length(res),
        2L
    )
    asub <- c(1,4)
    res <- ExpList[, , asub]
    expect_identical(
        length(res),
        2L
    )
    asub <- c(TRUE, FALSE)
    res <- ExpList[, , asub]
    expect_identical(
        length(res),
        2L
    )
})
