context("Colname consistency after subsetting")

test_that("when columns in different order, subset ops will rearrange them", {
    sampMap <- sampleMap(mae)
    ## Reverse the order of the samples in the first assay
    newSampMap <- rbind(sampMap[4:1, ], sampMap[5:nrow(sampMap), ])

    sampleMap(mae) <- newSampMap

    # Rearrange ExperimentList colnames using sampleMap "assay" column
    mae <- mae[, rownames(colData(mae)), ]

    expect_true(identical(sampleMap(mae)[1:4, "colname"], colnames(mae)[[1]]))
})
