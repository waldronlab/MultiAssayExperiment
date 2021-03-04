## Run the example ExperimentList
example("ExperimentList")

## Create sample maps for each experiment
exprmap <- data.frame(
    primary = c("Jack", "Jill", "Barbara", "Bob"),
    colname = c("array1", "array2", "array3", "array4"),
    stringsAsFactors = FALSE)

methylmap <- data.frame(
    primary = c("Jack", "Jack", "Jill", "Barbara", "Bob"),
    colname = c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5"),
    stringsAsFactors = FALSE)

rnamap <- data.frame(
    primary = c("Jack", "Jill", "Bob", "Barbara"),
    colname = c("samparray1", "samparray2", "samparray3", "samparray4"),
    stringsAsFactors = FALSE)

gistmap <- data.frame(
    primary = c("Jack", "Bob", "Jill"),
    colname = c("samp0", "samp1", "samp2"),
    stringsAsFactors = FALSE)

## Combine as a named list and convert to a DataFrame
maplist <- list(Affy = exprmap, Methyl450k = methylmap,
    RNASeqGene = rnamap, GISTIC = gistmap)

## Create a sampleMap
sampMap <- listToMap(maplist)
## Create an example phenotype data
colDat <- data.frame(sex = c("M", "F", "M", "F"), age = 38:41,
    row.names = c("Jack", "Jill", "Bob", "Barbara"))

## Create a MultiAssayExperiment instance
mae <- MultiAssayExperiment(experiments = ExpList, colData = colDat,
    sampleMap = sampMap)
