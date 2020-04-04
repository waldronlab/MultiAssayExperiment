context("check MultiAssayExperiment exportClass")

test_that("exportClass on a MultiAssayExperiment works", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]
    metadata(miniACC)[["exampleObj"]] <-
        GRanges(c(a = "chr2:1-10:-", b = "chr2:2-10:+", c = "chr2:3-10:+"))

    .write.table <- function(...) {
        return(..2)
    }

    with_mock(write.table = .write.table, {
        filenames <- basename(
            exportClass(miniACC, dir = tempdir(), fmt = "csv", ext = ".csv",
                verbose = FALSE)
        )
        expect_match(
            filenames, "miniACC_.*"
        )
    })

    with_mock(write.table = .write.table, {
        filenames <- basename(
            .sortMetadata(
                miniACC, "miniACC", dir = tempdir(), fmt = "csv", ext = ".csv"
            )
        )
        expect_match(
            filenames, "miniACC_META_.*"
        )
    })
})

test_that("constructors work", {
    expect_true(validObject(MultiAssayExperiment()))
    expect_true(validObject(MultiAssayExperiment(ExperimentList())))
    expect_true(validObject(MatchedAssayExperiment()))
    expect_true(validObject(MatchedAssayExperiment(ExperimentList())))
    ## remove replicate column
    expect_true(validObject(
        MatchedAssayExperiment(
            mae[, list(Affy = 1:4, Methyl450k = 2:5,
                    RNASeqGene = 1:4, GISTIC = 1:3), ]
        )
    ))
    expect_error(
        MatchedAssayExperiment(ExpList, colDat, sampMap)
    )
    mm <- mergeReplicates(mae)
    conv <- as(mm, "MatchedAssayExperiment")
    expect_true(
        is(conv, "MatchedAssayExperiment")
    )
})

test_that("replace methods are using rebliss and replace", {
    mae0 <- mae
    sampleMap(mae0) <- DataFrame()
    expect_true(isEmpty(mae0))
    expect_true(validObject(mae0))

    mae0 <- mae
    sampleMap(mae0) <- DataFrame(assay="testAssay", primary="testPrimary",
        colname="testColname")
    expect_true(isEmpty(mae0))
    expect_true(validObject(mae0))


    mae0 <- mae
    experiments(mae0) <- ExperimentList()
    expect_true(isEmpty(mae0))
    expect_true(validObject(mae0))

    mae0 <- mae
    colData(mae0) <- DataFrame()
    # check cols are zero
    expect_true(
        all(vapply(experiments(mae0),
            function(expo) dim(expo)[[2]] == 0L,
            logical(1)
        ))
    )
    expect_true(validObject(mae0))
})
