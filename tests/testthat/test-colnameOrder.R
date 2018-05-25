context("Colname consistency after subsetting")

test_that("when columns in different order, subset ops will rearrange them", {
    example("MultiAssayExperiment")
    sampMap <- sampleMap(myMultiAssayExperiment)
    ## Reverse the order of the samples in the first assay
    newSampMap <- rbind(sampMap[4:1, ], sampMap[5:nrow(sampMap), ])

    sampleMap(myMultiAssayExperiment) <- newSampMap

    # Rearrange ExperimentList colnames using sampleMap "assay" column
    myMultiAssayExperiment <-
        myMultiAssayExperiment[, rownames(colData(myMultiAssayExperiment)), ]

    expect_true(identical(sampleMap(myMultiAssayExperiment)[1:4, "colname"],
        colnames(myMultiAssayExperiment)[[1]]))
})
