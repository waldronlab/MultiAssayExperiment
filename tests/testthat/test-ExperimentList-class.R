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

test_that("Metadata is kept in ExperimentList when replacing", {
    example("MultiAssayExperiment")
    ## add metadata columns and metadata
    mcoldf <- DataFrame(AssayNumber=seq_len(length(ExpList)),
        row.names = names(ExpList))
    mcols(ExpList) <- mcoldf

    metalist <- list(Shiny = "Blue Jeans", Old = "Metadata")
    metadata(ExpList) <- metalist

    experiments(mae) <- ExpList

    expect_identical(mcols(experiments(mae)), mcoldf)
    expect_identical(metadata(experiments(mae)), metalist)
})
