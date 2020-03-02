context("check MultiAssayExperiment exportClass")

test_that("exportClass on a MultiAssayExperiment works", {
    env <- new.env(parent = emptyenv())
    data("miniACC",  envir = env)
    miniACC <- env[["miniACC"]]
    metadata(miniACC)[["exampleObj"]] <-
        GRanges(c(a = "chr2:1-10:-", b = "chr2:2-10:+", c = "chr2:3-10:+"))

    .write.table <- function(...) {
        message(..2)
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
