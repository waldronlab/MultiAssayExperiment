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
    expect_true(
        all.equal(
        vapply(replicated(obs2), function(x) unique(lengths(x)), integer(1L)),
        lengths(Filter(length, colnames(obs2)))
        )
    )
    expect_identical(
        showReplicated(obs2),
        list(
            m = CharacterList(structure(list(), .Names = character(0L))),
            m3 = CharacterList(structure(list(), .Names = character(0L))),
            m4 = CharacterList(a = c("A", "B"))
        )
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


test_that("merging replicates orders columns of each experiment the same", {

    check_mergeReps_colorder <- function(mae) {
        maeredux <- mergeReplicates(intersectColumns(mae))
        maeassay <- assays(maeredux)
        matchednames <- lapply(maeassay, function(x) {
            sampleMap(maeredux)[
                match(colnames(x), sampleMap(maeredux)$colname),
                "primary"
            ]
        })
        for (i in seq_along(matchednames)[-1])
            expect_identical(matchednames[[1]], matchednames[[i]])
    }

    mae0 <- mae
    check_mergeReps_colorder(mae0)

    mae0[[1]] <- mae[[1]][, 4:1]
    check_mergeReps_colorder(mae0)
    check_mergeReps_colorder(mae0[, , 1:2])
    check_mergeReps_colorder(mae0[, , 1:3])
})


test_that("getWithColData works", {
    Y1 <- matrix(rnorm(200), ncol=10)
    colnames(Y1) <- c("A.1", "A.2",
        "B.1", "B.2", "B.3",
        "C.1", "C.2",
        "D.1",
        "E.1", "E.2")

    # Making up some sample-level metadata, for some variety.
    df <- DataFrame(RIN=runif(10), sizeFactor=runif(10))
    rownames(df) <- colnames(Y1)
    se1 <- SummarizedExperiment(list(counts=Y1), colData=df)

    ### EXPERIMENT 2 ###
    Y2 <- matrix(rnorm(200), ncol=5)
    colnames(Y2) <- c("B", "C", "D", "E", "F")

    # Making up some more sample-level metadata.
    df <- DataFrame(FrIP=runif(5), NumPeaks=sample(1000, 5))
    rownames(df) <- colnames(Y2)
    se2 <- SummarizedExperiment(list(counts=Y2), colData=df)

    sampleMap <- rbind(
        DataFrame(assay="rnaseq",
            primary=sub("\\..*", "", colnames(se1)),
            colname=colnames(se1)
        ),
        DataFrame(assay="chipseq",
            primary=sub("\\..*", "", colnames(se2)),
            colname=colnames(se2)
        )
    )

    colData <- DataFrame(
        row.names=LETTERS[1:6],
        Sex=sample(c("M", "F"), 6, replace=TRUE),
        Age=sample(10:50, 6)
    )

    MAE <- MultiAssayExperiment(list(rnaseq=se1, chipseq=se2),
        sampleMap=sampleMap, colData=colData)

    expect_warning(
        getWithColData(MAE, 1L, "append"),
        "^Duplicating"
    )

    expect_warning(
        getWithColData(MAE, 1L, "replace"),
        "^Duplicating"
    )

    cData <- colData(getWithColData(MAE, 1L, "append"))
    cDataMatch <- cData[, names(colData), drop = FALSE]
    matchedData <- colData[
        match(mapToList(sampleMap)[[1]][["primary"]], rownames(colData)), ,
        drop = FALSE]
    testres <- unlist(
        Map(function(x, y) { identical(x, y) }, x = cDataMatch, y = matchedData)
    )
    expect_true( all(testres) )

    expect_identical(
        names(cData),
        c(names(colData(MAE[[1]])), names(colData))
    )
    expect_identical(
        colData(se1),
        cData[, names(colData(se1))]
    )
    testres <- unlist(
        Map(function(x, y) { identical(x, y) },
            x = colData(se1), y = cData[, names(colData(se1))])
    )
    expect_true( all(testres) )

    cData <- colData(getWithColData(MAE, 1L, "replace"))
    cDataMatch <- cData[, names(colData), drop = FALSE]
    matchedData <- colData[
        match(mapToList(sampleMap)[[1L]][["primary"]], rownames(colData)), ,
        drop = FALSE]
    expect_identical(cDataMatch, matchedData)

    ## Empty MAE colData return assay colData for append
    MAE0 <- MAE
    colData(MAE0) <- DataFrame(row.names = rownames(colData(MAE0)))
    cData <- colData(getWithColData(MAE0, 1L, "append"))
    expect_identical(
        cData,
        colData(se1)
    )

    ## Empty MAE colData return assay colData for replace
    MAE0 <- MAE
    colData(MAE0) <- DataFrame(row.names = rownames(colData(MAE0)))
    cData <- colData(getWithColData(MAE0, 1L, "replace"))

    expect_identical(
        cData,
        DataFrame(row.names = mapToList(sampleMap(MAE0))[[1]][["primary"]])
    )
    expect_true( isEmpty(cData) )


    ## Shuffled rows in MAE colData for append
    MAE0 <- MAE
    colData(MAE0) <- colData(MAE0)[sample(nrow(colData(MAE0))), ]
    cData <- colData(getWithColData(MAE0, 1L, "append"))
    cDataMatch <- cData[, names(colData(MAE0)), drop = FALSE]
    matchedData <- colData(MAE0)[
        match(
            mapToList(sampleMap)[[1]][["primary"]], rownames(colData(MAE0))
        ), , drop = FALSE]
    matchedData <- matchedData[order(rownames(matchedData)), ]
    testres <- unlist(
        Map(function(x, y) { identical(x, y) }, x = cDataMatch, y = matchedData)
    )
    expect_true( all(testres) )

    ## Shuffled rows in MAE colData for replace
    MAE0 <- MAE
    colData(MAE0) <- colData(MAE0)[sample(nrow(colData(MAE0))), ]
    cData <- colData(getWithColData(MAE0, 1L, "replace"))
    matchedData <- colData(MAE0)[
        match(
            mapToList(sampleMap)[[1]][["primary"]], rownames(colData(MAE0))
        ), , drop = FALSE]
    expect_identical(cData, matchedData)

    ## Empty assay colData return MAE colData
    MAE0 <- MAE
    eCol <- DataFrame(row.names = rownames(colData(MAE0[["rnaseq"]])))
    colData(MAE0[["rnaseq"]]) <- eCol
    cData <- colData(getWithColData(MAE0, 1L, "replace"))
    matchedData <- colData(MAE0)[
        match(mapToList(sampleMap(MAE0))[["rnaseq"]][["primary"]],
            rownames(colData)), ,
        drop = FALSE]
    expect_identical(cData, matchedData)


    MAE0 <- MAE
    eCol <- DataFrame(row.names = rownames(colData(MAE0[["rnaseq"]])))
    colData(MAE0[["rnaseq"]]) <- eCol
    cData <- colData(getWithColData(MAE0, 1L, "append"))
    matchedData <- colData(MAE0)[
        match(mapToList(sampleMap(MAE0))[["rnaseq"]][["primary"]],
            rownames(colData)), ,
        drop = FALSE]
    expect_identical(cData, matchedData)

    ## Test for row order mismatch bug (Issue: James Bonaffini)
    ## When experiment colData row order differs from sampleMap order,
    ## appended MAE colData should still be correctly matched
    Y3 <- matrix(rnorm(40), ncol=4)
    colnames(Y3) <- c("S001", "S002", "S003", "S004")

    # Experiment colData with specific row order and a common column
    exp_df <- DataFrame(PATNUM=c("P3", "P1", "P4", "P2"),
                        SAMPLE_TYPE=c("tumor", "normal", "tumor", "normal"),
                        row.names = colnames(Y3))
    se3 <- SummarizedExperiment(list(counts=Y3), colData=exp_df)

    # MAE colData with patient-level metadata
    mae_colData <- DataFrame(PATNUM=c("P1", "P2", "P3", "P4"),
                             AGE=c(45, 67, 52, 38),
                             TREATMENT=c("A", "B", "C", "D"),
                             row.names=c("P1", "P2", "P3", "P4"))

    # SampleMap with different order than experiment colData
    smap <- DataFrame(assay=factor("exp3"),
                      primary=c("P1", "P2", "P3", "P4"),
                      colname=c("S002", "S004", "S001", "S003"))

    MAE_test <- MultiAssayExperiment(list(exp3=se3),
                                     sampleMap=smap,
                                     colData=mae_colData)

    result <- getWithColData(MAE_test, 1L, "append")
    resultData <- colData(result)

    # Verify each sample has correct patient's metadata
    expect_identical(resultData["S001", "PATNUM"], "P3")
    expect_identical(resultData["S001", "AGE"], 52)
    expect_identical(resultData["S001", "TREATMENT"], "C")

    expect_identical(resultData["S002", "PATNUM"], "P1")
    expect_identical(resultData["S002", "AGE"], 45)
    expect_identical(resultData["S002", "TREATMENT"], "A")

    expect_identical(resultData["S003", "PATNUM"], "P4")
    expect_identical(resultData["S003", "AGE"], 38)
    expect_identical(resultData["S003", "TREATMENT"], "D")

    expect_identical(resultData["S004", "PATNUM"], "P2")
    expect_identical(resultData["S004", "AGE"], 67)
    expect_identical(resultData["S004", "TREATMENT"], "B")

    # Verify MAE-only columns match correctly for all rows
    mae_only_cols <- c("AGE", "TREATMENT")
    for (i in seq_len(nrow(resultData))) {
        sample_id <- rownames(resultData)[i]
        patient_id <- resultData[i, "PATNUM"]
        expected_data <- mae_colData[patient_id, mae_only_cols, drop = FALSE]
        actual_data <- resultData[i, mae_only_cols, drop = FALSE]
        testres <- unlist(
            Map(function(x, y) { identical(x, y) }, x=actual_data, y=expected_data)
        )
        expect_true(all(testres))
    }
})

test_that("renaming helpers work", {
    ## renamePrimary
    newnames <- c("Red", "Green", "Blue", "Yellow")
    mae0 <- renamePrimary(mae, newnames)
    expect_identical(
        rownames(colData(mae0)),
        newnames
    )
    expect_true(
        all(newnames %in% sampleMap(mae0)[["primary"]])
    )

    ## renameColname
    newaffynames <- paste0("Affy", 1:4)
    mae0 <- renameColname(mae, "Affy", newaffynames)
    expect_identical(
        colnames(mae0)[["Affy"]],
        newaffynames
    )
})


test_that("longForm,ANY-method works", {
    denv <- new.env(parent = emptyenv())
    data("miniACC", package="MultiAssayExperiment", envir = denv)
    miniACC <- denv[["miniACC"]]

    expect_silent(
        longdf <- longForm(miniACC)
    )
    expect_true(
        all(
            as.character(longdf[["rowname"]]) %in%
                Reduce(union, rownames(miniACC))
        )
    )

    m <- matrix(1:12, nrow = 4)
    colnames(m) <- LETTERS[1:3]
    se <- SummarizedExperiment(m)
    rownames(se) <- c(1, 2, 4, 5)

    mae <- MultiAssayExperiment(list(x = se))
    longform <- longForm(mae)
    expect_true(
        is.character(longform[["rowname"]])
    )
    expect_identical(
        typeof(m),
        typeof(longform[["value"]])
    )
})
