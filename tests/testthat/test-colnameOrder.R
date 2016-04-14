context("Colname consistency after subsetting")

example("MultiAssayExperiment")

## Reverse the order of the samples in the first assay
newSampMap <- rbind(sampleMap(myMultiAssayExperiment)[4:1, ],
                    sampleMap(myMultiAssayExperiment)[
                      5:nrow(sampleMap(myMultiAssayExperiment)),
                      ])
sampleMap(myMultiAssayExperiment) <- newSampMap

# Rearrange Elist colnames using sampleMap "assay" column
myMultiAssayExperiment <-
  myMultiAssayExperiment[, rownames(pData(myMultiAssayExperiment)), ]


test_that("when columns in different order, subset ops will rearrange them", {
  expect_true(identical(sampleMap(myMultiAssayExperiment)[1:4, "assay"],
                        colnames(myMultiAssayExperiment)[[1]]))
})
