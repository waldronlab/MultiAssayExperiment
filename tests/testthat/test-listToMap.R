context("listToMap converts a DataFrame list to sampleMap")

test_that("listToMap works with data.frame", {
    samp <- data.frame(
        assay = rep(letters[1:3], each = 3),
        primary = rep(paste0(letters[1:3], 1:3), each = 3),
        colname = as.character(1:9)
    )

    ssamp <- split(samp[-1], samp[["assay"]])
    Dsamp <- as(samp, "DataFrame")
    Dsamp[["assay"]] <- as.factor(Dsamp[["assay"]])

    identical(Dsamp, listToMap(ssamp))
})

test_that("listToMap works with DataFrame", {
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

test_that("listToMap works with SplitDataFrameList", {
    samp <- data.frame(
        assay = rep(letters[1:3], each = 3),
        primary = rep(paste0(letters[1:3], 1:3), each = 3),
        colname = as.character(1:9)
    )

    ssamp <- split(samp[-1], samp[["assay"]])
    SDFLsamp <- as(ssamp, "SplitDataFrameList")

    DFsamp <- as(samp, "DataFrame")
    DFsamp[["assay"]] <- as.factor(DFsamp[["assay"]])

    expect_identical(DFsamp, listToMap(SDFLsamp))
})

test_that("listToMap fill works with empty elements", {
    samps <- paste0("samp", 1:2)
    listMap <- list(
        a1 = DataFrame(primary = samps, colname = samps),
        a2 = DataFrame(primary = character(0L), colname = character(0L))
    )
    nmap <- listToMap(listMap, fill = FALSE)
    expect_identical(
        nmap,
        DataFrame(
            assay = factor(c("a1", "a1"), levels = c("a1", "a2")),
            primary = samps,
            colname = samps
        )
    )
    nmap <- listToMap(listMap, fill = TRUE)
    expect_identical(
        nmap,
        DataFrame(
            assay = factor(c("a1", "a1", "a2"), levels = c("a1", "a2")),
            primary = c(samps, NA_character_),
            colname = c(samps, NA_character_)
        )
    )
})

test_that("listToMap uses list names as factor levels", {
    listMap <- list(
        assay1 =
            DataFrame(primary = paste0("A", 1:4), colname = paste0("A", 1:4)),
        assay2 =
            DataFrame(primary = character(0L), colname = character(0L))
    )
    result <- listToMap(listMap, fill = FALSE)
    expect_identical(
        levels(result[["assay"]]),
        c("assay1", "assay2")
    )
})
