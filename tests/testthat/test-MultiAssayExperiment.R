test_that("MultiAssayExperiment constructor works", {

    expect_true(validObject(MultiAssayExperiment()))

    expect_error(MultiAssayExperiment(list(matrix(0, 0, 0))),
                 label="unnamed ExperimentList")
    
    obs <- MultiAssayExperiment(list(m=matrix(0, 0, 0)))
    expect_true(validObject(obs))
    expect_identical(names(experiments(obs)), "m")
    
    pData = S4Vectors::DataFrame(row.names=letters[1:4])
    expect_true(validObject(MultiAssayExperiment(pData=pData)))

    sampleMap = S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs <- MultiAssayExperiment(sampleMap=sampleMap)
    expect_true(validObject(obs))
    expect_identical(nrow(sampleMap(obs)), 0L)
    expect_identical(length(experiments(obs)), 0L)
    expect_identical(dim(pData(obs)), c(0L, 0L))

  
})


test_that("MultiAssayExperiment .harmonize construction helper works", {

    
    m = matrix(0, 2, 2, dimnames=list(letters[1:2], letters[1:2]))
    experiments = ExperimentList(list(m=m))
    sampleMap = S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary="a", colname="a")
    obs = MultiAssayExperiment(experiments, sampleMap=sampleMap)
    # remove unused assays in experiments not in sampleMap 
    expect_identical(colnames(obs)[[1]], "a")
    
    sampleMap = S4Vectors::DataFrame(assay=factor("m", levels="m"),
        primary=letters, colname=letters)
    obs = MultiAssayExperiment(experiments, sampleMap=sampleMap)
    # removed unused experiments in sampleMap
    expect_identical(colnames(experiments)[[1]], sampleMap(obs)[["colname"]])
    # pData initialized when not specified
    expect_identical(colnames(experiments)[[1]], rownames(pData(obs)))
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

    # pData subset by sampleData primary 
    pData = S4Vectors::DataFrame(matrix(0, 4, 4, dimnames=list(letters[1:4], letters[1:4])))
    obs = MultiAssayExperiment(experiments, sampleMap=sampleMap, pData=pData)
    expect_identical(rownames(pData(obs)), sampleMap(obs)[["primary"]])

    # combo
    sampleMap = S4Vectors::DataFrame(
        assay=factor(c("m","stack"), levels=c("m","stack")),
        primary=c("a","n"),
        colname=c("a","n"))
    m=matrix(0, 3, 3, dimnames=list(letters[1:3], letters[1:3]))
    experiments=ExperimentList(list(m=m))
    obs = MultiAssayExperiment(experiments=experiments, sampleMap=sampleMap)
    
}
