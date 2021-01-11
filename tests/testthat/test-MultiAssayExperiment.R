test_that("MultiAssayExperiment constructor works", {

    obs <- MultiAssayExperiment()
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), character())
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(colData(obs)), c(0L, 0L))

    expect_error(MultiAssayExperiment(list(matrix(0, 0, 0))),
                 label="unnamed ExperimentList")

    obs <- MultiAssayExperiment(list(m=matrix(0, 0, 0)))
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), "m")
    expect_identical(length(experiments(obs)), 1L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(colData(obs)), c(0L, 0L))

    colData <- S4Vectors::DataFrame(row.names=letters[1:4])
    obs <- MultiAssayExperiment(colData=colData)
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), character())
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(colData(obs)), c(0L, 0L))

    sampleMap <- S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs <- MultiAssayExperiment(sampleMap=sampleMap)
    expect_true(validObject(obs))
    expect_identical(nrow(sampleMap(obs)), 0L)
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(colData(obs)), c(0L, 0L))

    # test multiple experiment
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 1, 1, dimnames=list("d", "d"))
    obs <- MultiAssayExperiment(list(m=m, m2=m2, m3=m3))
    expect_true(validObject(obs))
    expect_identical(length(experiments(obs)), 3L)
    expect_identical(levels(sampleMap(obs)[["assay"]]), c("m", "m2", "m3"))
    expect_identical(dim(sampleMap(obs)), c(4L, 3L))
    expect_identical(rownames(colData(obs)), sampleMap(obs)[["primary"]])

})


test_that("MultiAssayExperiment .harmonize construction helper works", {

    # remove unused assays in experiments not in sampleMap
    m <- matrix(0, 2, 2, dimnames=list(letters[1:2], letters[1:2]))
    experiments <- ExperimentList(list(m=m))
    sampleMap <- S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary="a", colname="a")
    obs <- MultiAssayExperiment(experiments, sampleMap=sampleMap)
    expect_identical(colnames(obs)[[1]], "a")
    expect_identical(sampleMap(obs), sampleMap)
    expect_identical(row.names(colData(obs)), "a")

    # removed unused assays in sampleMap
    sampleMap <- S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs <- MultiAssayExperiment(experiments, sampleMap=sampleMap)
    expect_identical(colnames(experiments)[[1]], sampleMap(obs)[["colname"]])
    expect_identical(experiments(obs), experiments)
    expect_identical(colnames(experiments)[[1]], rownames(colData(obs)))
    expect_identical(rownames(colData(obs)), sampleMap(obs)[["primary"]])

    # colData subset by sampleData primary
    colData <- S4Vectors::DataFrame(matrix(0, 4, 4, dimnames=list(letters[1:4],
                                                                letters[1:4])))
    obs <- MultiAssayExperiment(experiments, sampleMap=sampleMap, colData=colData)
    expect_identical(rownames(colData(obs)), sampleMap(obs)[["primary"]])

    # combo subset
    sampleMap <- S4Vectors::DataFrame(
        assay=factor(c("m", "stack"), levels=c("m", "stack")),
        primary=c("a", "n"),
        colname=c("a", "n"))
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    colData <- S4Vectors::DataFrame(matrix(0, 4, 4, dimnames=list(letters[1:4],
                                                                letters[1:4])))
    experiments <- ExperimentList(list(m=m))
    obs <- MultiAssayExperiment(experiments=experiments,
                                sampleMap=sampleMap, colData=colData)
    expect_identical(colnames(obs)[[1]], "a")
    expect_identical(dim(sampleMap(obs)), c(1L, 3L))
    expect_identical(rownames(colData(obs)), sampleMap(obs)[["primary"]])

    # test multiple experiment, combo subsetting
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 2, 2, dimnames=list(letters[4:5], letters[4:5]))
    experiments <- ExperimentList(list(m=m, m2=m2, m3=m3))
    sampleMap <- S4Vectors::DataFrame(
        assay=factor(c("m", "stack", "m3"), levels=c("m", "stack", "m3")),
        primary=c("a", "n", "d"),
        colname=c("a", "n", "d"))
    colData <- S4Vectors::DataFrame(matrix(0, 7, 7,
                                         dimnames=list(letters[1:7],
                                                       letters[1:7])))
    obs <- MultiAssayExperiment(experiments=experiments,
                                sampleMap=sampleMap,
                                colData=colData)
    expect_true(validObject(obs))
    expect_identical(length(experiments(obs)), 2L)
    expect_identical(dim(experiments(obs)[[1]]), c(3L, 1L))
    expect_identical(dim(experiments(obs)[[2]]), c(2L, 1L))
    expect_identical(names(experiments(obs)), c("m", "m3"))
    expect_identical(dim(sampleMap(obs)), c(2L, 3L))
    expect_identical(levels(sampleMap(obs)[["assay"]]), c("m", "m3"))
    expect_identical(rownames(colData(obs)), sampleMap(obs)[["primary"]])
    expect_identical(unlist(colnames(experiments(obs)), use.names=FALSE),
                     row.names(colData(obs)))

})

test_that("MultiAssayExperiment replacements work", {
    pDF <- DataFrame(a = 1:4, b = letters[1:4])
    obs <- MultiAssayExperiment()
    expect_error(colData(obs) <- pDF)
})

test_that("MultiAssayExperiment name replacements work", {
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 1, 1, dimnames=list("d", "d"))
    obs <- MultiAssayExperiment(list(m=m, m2=m2, m3=m3))
    newnames <- c("exp1", "exp2", "exp3")
    names(obs) <- newnames
    expect_identical(names(obs), newnames)
    expect_identical(levels(sampleMap(obs)[["assay"]]), newnames)
})

test_that("ExperimentList metadata and mcols are preserved", {
    mcoldf <- DataFrame(AssayNumber=seq_len(length(ExpList)),
        row.names = names(ExpList))
    mcols(ExpList) <- mcoldf

    metalist <- list(Shiny = "Blue Jeans", Old = "Metadata")
    metadata(ExpList) <- metalist

    mae0 <- MultiAssayExperiment(ExpList)

    expect_identical(mcols(experiments(mae0)), mcoldf)
    expect_identical(metadata(experiments(mae0)), metalist)
})
