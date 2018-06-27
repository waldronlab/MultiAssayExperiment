## Load example MultiAssayExperiment
example(MultiAssayExperiment)

## Access the sampleMap
sampleMap(mae)

## Replacement method for a MultiAssayExperiment sampleMap
sampleMap(mae) <- S4Vectors::DataFrame()

## Access the ExperimentList
experiments(mae)

## Replace with an empty ExperimentList
experiments(mae) <- ExperimentList()

## Access the metadata
metadata(mae)

## Replace metadata with a list
metadata(mae) <- list(runDate =
    format(Sys.time(), "%B %d, %Y"))

## Access the colData
colData(mae)

## Access a column in colData
mae$age

## Replace a column in colData
mae$age <- mae$age + 1
