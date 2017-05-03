## Load example MultiAssayExperiment
example(MultiAssayExperiment)

## Access the sampleMap
sampleMap(myMultiAssayExperiment)

## Replacement method for a MultiAssayExperiment sampleMap
sampleMap(myMultiAssayExperiment) <- S4Vectors::DataFrame()

## Access the ExperimentList
experiments(myMultiAssayExperiment)

## Replace with an empty ExperimentList
experiments(myMultiAssayExperiment) <- ExperimentList()

## Access the metadata
metadata(myMultiAssayExperiment)

## Replace metadata with a list
metadata(myMultiAssayExperiment) <- list(runDate =
                                             format(Sys.time(), "%B %d, %Y"))

## Access a column in colData
myMultiAssayExperiment$age

## Replace a column in colData
myMultiAssayExperiment$age <- myMultiAssayExperiment$age + 1
