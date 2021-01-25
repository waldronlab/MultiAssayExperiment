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

    mat <- matrix(rpois(10000, 5), nrow=10)
    rownames(mat) <- paste0("CD", seq_len(nrow(mat))) # made-up markers
    colnames(mat) <- paste0("BARCODE_", seq_len(ncol(mat))) # made-up barcodes

    # barcodes come from a range of primaries
    primaries <- rep(LETTERS[1:10], each=100)

    # test constructor function works with only sampleMap and colData info
    mae <- MultiAssayExperiment(
        experiments = list(markers = mat),
        sampleMap = data.frame(
            assay = "markers", primary = primaries, colname = colnames(mat)
        ),
        colData=DataFrame(row.names=LETTERS[1:10])
    )
    expect_true(validObject(mae))
    expect_equal(unique(sampleMap(mae)[["primary"]]), unique(primaries))

    mae <- MultiAssayExperiment(
        experiments = list(markers = mat),
        sampleMap = data.frame(
            assay = "markers", primary = primaries, colname = colnames(mat)
        ),
    )
    expect_true(validObject(mae))
    expect_equal(sampleMap(mae)[["primary"]], primaries)


    exprss1 <- matrix(rnorm(16), ncol = 4,
        dimnames = list(sprintf("ENST00000%i", sample(288754:290000, 4)),
            c("Jack", "Jill", "Bob", "Bobby"))
    )
    exprss2 <- matrix(rnorm(12), ncol = 3,
        dimnames = list(sprintf("ENST00000%i", sample(288754:290000, 4)),
            c("Jack", "Jane", "Bob"))
    )
    expl <- list("methyl 2k"  = exprss1, "methyl 3k" = exprss2)
    ptdata <- data.frame(
        sex=c("M", "F", "M", "F"),
        age=38:41,
        row.names=c("Jack", "Jill", "Bob", "Barbara")
    )

    mae <- MultiAssayExperiment(experiments=expl, colData=ptdata)
    expect_true(validObject(mae))
    expect_true(!isEmpty(mae))
    expect_true(length(mae) == 2L)
    expect_equivalent(vapply(experiments(mae), nrow, integer(1L)), c(4L, 4L))
    expect_equivalent(vapply(experiments(mae), ncol, integer(1L)), c(3L, 2L))
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

    mae0 <- mae
    expect_error(
        colData(mae0) <- DataFrame(rownames = "Blue")
    )

    cc <- colnames(mae)
    cc[[1]] <- toupper(cc[[1]])
    cc[[3]] <- toupper(cc[[3]])

    mae0 <- mae
    colnames(mae0) <- cc
    expect_identical(colnames(mae0), cc)

    mae0 <- mae
    cc <- as.list(cc)
    colnames(mae0) <- cc
    expect_identical(
        colnames(mae0), as(cc, "CompressedCharacterList")
    )
})
