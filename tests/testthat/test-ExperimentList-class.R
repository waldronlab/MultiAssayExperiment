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

test_that("ExperimentList validity function catches non-compatible classes", {
    # list
    testObj <- list(A = 1:3, B = letters[1:3])
    expect_error(ExperimentList(assay1 = testObj))
    # List
    testObj <- List(A = 1:3, B = letters[1:3])
    expect_error(ExperimentList(assay1 = testObj))
    testObj <- list(A = matrix(1, 1, dimnames = list("A", "A")))
    expect_true(validObject(ExperimentList(assay1 = testObj)))
    A <- array(1:24, 4:2, dimnames =
        list(letters[1:4], LETTERS[1:3], c("A", "B"))
    )
    B <- matrix(1:12, ncol=3, dimnames = list(letters[1:4], LETTERS[1:3]))
    se2 <- SummarizedExperiment(list(A=A, B=B))
    expect_true(validObject(ExperimentList(se2 = se2)))
    # top-level DataFrame / data.frame
    testObj <- DataFrame(A = 1:3, B = letters[1:3])
    testObj1 <- data.frame(A = 1:3, B = letters[1:3])
    expect_warning(ExperimentList(assay1 = list(a = testObj)))
    expect_warning(ExperimentList(assay1 = list(a = testObj, b = testObj1)))
    expect_warning(ExperimentList(assay1 = testObj))
    ## GRangesList
    expect_error(
        ExperimentList(list(GR = GRangesList())),
        "GRangesList"
    )
    ## vector
    expect_error(
        ExperimentList(list(int = 1:100L)),
        "integer"
    )
})
