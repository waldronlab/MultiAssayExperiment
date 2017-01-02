test_that("MultiAssayExperiment constructor works", {

    obs <- MultiAssayExperiment()
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), character())
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(pData(obs)), c(0L, 0L))

    expect_error(MultiAssayExperiment(list(matrix(0, 0, 0))),
                 label="unnamed ExperimentList")

    obs <- MultiAssayExperiment(list(m=matrix(0, 0, 0)))
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), "m")
    expect_identical(length(experiments(obs)), 1L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(pData(obs)), c(0L, 0L))

    pData <- S4Vectors::DataFrame(row.names=letters[1:4])
    obs <- MultiAssayExperiment(pData=pData)
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), character())
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(sampleMap(obs)), c(0L, 3L))
    expect_identical(dim(pData(obs)), c(0L, 0L))

    sampleMap <- S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs <- MultiAssayExperiment(sampleMap=sampleMap)
    expect_true(validObject(obs))
    expect_identical(nrow(sampleMap(obs)), 0L)
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(pData(obs)), c(0L, 0L))

    # test multiple experiment
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 1, 1, dimnames=list("d", "d"))
    obs <- MultiAssayExperiment(list(m=m, m2=m2, m3=m3))
    expect_true(validObject(obs))
    expect_identical(length(experiments(obs)), 3L)
    expect_identical(levels(sampleMap(obs)[["assay"]]), c("m", "m2", "m3"))
    expect_identical(dim(sampleMap(obs)), c(4L, 3L))
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

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
    expect_identical(row.names(pData(obs)), "a")

    # removed unused assays in sampleMap
    sampleMap <- S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs <- MultiAssayExperiment(experiments, sampleMap=sampleMap)
    expect_identical(colnames(experiments)[[1]], sampleMap(obs)[["colname"]])
    expect_identical(experiments(obs), experiments)
    expect_identical(colnames(experiments)[[1]], rownames(pData(obs)))
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

    # pData subset by sampleData primary
    pData <- S4Vectors::DataFrame(matrix(0, 4, 4, dimnames=list(letters[1:4],
                                                                letters[1:4])))
    obs <- MultiAssayExperiment(experiments, sampleMap=sampleMap, pData=pData)
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

    # combo subset
    sampleMap <- S4Vectors::DataFrame(
        assay=factor(c("m", "stack"), levels=c("m", "stack")),
        primary=c("a", "n"),
        colname=c("a", "n"))
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    pData <- S4Vectors::DataFrame(matrix(0, 4, 4, dimnames=list(letters[1:4],
                                                                letters[1:4])))
    experiments <- ExperimentList(list(m=m))
    obs <- MultiAssayExperiment(experiments=experiments,
                                sampleMap=sampleMap, pData=pData)
    expect_identical(colnames(obs)[[1]], "a")
    expect_identical(dim(sampleMap(obs)), c(1L, 3L))
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

    # test multiple experiment, combo subsetting
    m <- matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    m2 <- matrix(0, 0, 0)
    m3 <- matrix(0, 2, 2, dimnames=list(letters[4:5], letters[4:5]))
    experiments <- ExperimentList(list(m=m, m2=m2, m3=m3))
    sampleMap <- S4Vectors::DataFrame(
        assay=factor(c("m", "stack", "m3"), levels=c("m", "stack", "m3")),
        primary=c("a", "n", "d"),
        colname=c("a", "n", "d"))
    pData <- S4Vectors::DataFrame(matrix(0, 7, 7,
                                         dimnames=list(letters[1:7],
                                                       letters[1:7])))
    obs <- MultiAssayExperiment(experiments=experiments,
                                sampleMap=sampleMap,
                                pData=pData)
    expect_true(validObject(obs))
    expect_identical(length(experiments(obs)), 2L)
    expect_identical(dim(experiments(obs)[[1]]), c(3L, 1L))
    expect_identical(dim(experiments(obs)[[2]]), c(2L, 1L))
    expect_identical(names(experiments(obs)), c("m", "m3"))
    expect_identical(dim(sampleMap(obs)), c(2L, 3L))
    expect_identical(levels(sampleMap(obs)[["assay"]]), c("m", "m3"))
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])
    expect_identical(unlist(colnames(experiments(obs)), use.names=FALSE),
                     row.names(pData(obs)))

})
