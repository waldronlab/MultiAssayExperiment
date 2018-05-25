context("sampleMap helpers")

test_that("listToMap works", {
    assayA <- DataFrame(primary = "pt1", colname = "samp1")
    assayB <- DataFrame(primary = "pt1", colname = "samp2")

    expect_identical(
        listToMap(list(A = assayA, B = assayB)),
        DataFrame(assay = factor(c("A", "B")),
            primary = c("pt1", "pt1"),
            colname = c("samp1", "samp2"))
    )

    expect_error(listToMap(list(DataFrame(letters[1:3], letters[3:1]))))

    newMap <- listToMap(
        list(
            assay1 = DataFrame(
                primary = paste0("pt", 1:3),
                colname = paste0("samp", 1:3)),
            assay2 = DataFrame(
                primary = paste0("pt", 1:3),
                colname = paste0("samp", 4:6))
        )
    )
    expect_identical(
        newMap,
        DataFrame(
            assay = factor(rep(c("assay1", "assay2"),
                times = c(3, 3))),
            primary = rep(paste0("pt", 1:3), 2),
            colname = paste0("samp", 1:6)
        )
    )
})
