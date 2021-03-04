example("MultiAssayExperiment")

## Create a sampleMap from a list using the listToMap function
sampMap <- listToMap(maplist)

## The inverse operation is also available
maplist <- mapToList(sampMap)
