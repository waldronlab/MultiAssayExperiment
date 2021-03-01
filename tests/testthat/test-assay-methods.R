context("assay methods work with ExperimentList when array present", {
    A <- array(1:24, 4:2, dimnames =
        list(letters[1:4], LETTERS[1:3], c("A", "B"))
    )
    se2 <- SummarizedExperiment(list(A=A))
    exp1 <- ExperimentList(A=se2)
    expect_true(validObject(exp1))
    # assay should return assay,SummarizedExperiment-method
    expect_true(is.array(assay(exp1)))

    A <- array(1:24, 4:2, dimnames =
        list(letters[1:4], LETTERS[1:3], c("A", "B"))
    )
    B <- matrix(1:12, ncol=3, dimnames = list(LETTERS[1:4], letters[1:3]))
    se2 <- SummarizedExperiment(list(A=A, B=B))
    exp1 <- ExperimentList(A=se2)
    expect_warning(assays(exp1))
    expect_true(is.array(assays(exp1)[[1]]))
})
