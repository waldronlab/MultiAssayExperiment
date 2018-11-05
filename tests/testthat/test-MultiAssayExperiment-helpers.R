test_that("MultiAssayExperiment anyReplicated returns same order", {
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 1, 1, dimnames=list("d", "d"))
    m4 <- matrix(0, 4, 4, dimnames=list(LETTERS[1:4], LETTERS[1:4]))
    obs <- MultiAssayExperiment(list(m=m, m2=m2, m3=m3))
    m4map <- DataFrame(
        assay = rep("m4", 4),
        primary = c("a", "a", "b", "d"),
        colname = LETTERS[1:4]
    )
    obs2 <- c(obs, m4 = m4, sampleMap = m4map)
    expect_true(
        length(replicated(obs2)[["m4"]]) == length(unique(m4map[["primary"]]))
    )
    expect_identical(
        unname(
        vapply(replicated(obs2), function(x) unique(lengths(x)), integer(1L))
        ),
        lengths(Filter(length, colnames(obs2)))
    )

    expect_identical(
        lapply(replicated(obs), names),
        lapply(mapToList(sampleMap(obs)), function(x) x[["primary"]])
    )

    ## test primary is identical to colnames
    cnames <- lapply(mapToList(sampleMap(obs)), `[[`, "colname")
    expect_identical(
        IRanges::CharacterList(cnames), Filter(length, colnames(obs))
    )
})
