## Load an example MultiAssayExperiment object
example("MultiAssayExperiment")
example("GRangesList")

## Find what ranges fit the criteria (see findOverlaps)
getHits(myMultiAssayExperiment, gr1)
