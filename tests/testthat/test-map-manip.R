context("sampleMap helpers")

test_that("listToMap works", {
    expect_identical(listToMap(list("B" = character(0L))),
    DataFrame(assay = factor("B"),
              primary = NA_character_,
              colname = NA_character_))
    expect_error(listToMap(list(DataFrame(letters[1:3], letters[3:1]))))
    newMap <- listToMap(list(
        assay1 = DataFrame(
            patient = paste0("pt", 1:3),
            sample = paste0("samp", 1:3)),
        assay2 = DataFrame(
            patient = paste0("pt", 1:3),
            sample = paste0("samp", 4:6))
    ))
    expect_identical(newMap, DataFrame(
        assay = factor(rep(c("assay1", "assay2"), times = c(3, 3))),
        primary = rep(paste0("pt", 1:3), 2),
        colname = paste0("samp", 1:6)))
})
