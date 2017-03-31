## Run the example ExperimentList
example("ExperimentList")

## Load example GRangesList object
example("RangedRaggedAssay")

## Add the RangedRaggedAssay to the list
ExpList <- c(ExpList, CNVgistic = myRRA)

## Create sample maps for each experiment
exprmap <- data.frame(
    primary = c("Jack", "Jill", "Barbara", "Bob"),
    colname = c("array1", "array2", "array3", "array4"),
    stringsAsFactors = FALSE)

methylmap <- data.frame(
    primary = c("Jack", "Jack", "Jill", "Barbara", "Bob"),
    colname = c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5"),
    stringsAsFactors = FALSE)

rangemap <- data.frame(
    primary = c("Jack", "Jill", "Jill"),
    colname = c("snparray1", "snparray2", "snparray3"),
    stringsAsFactors = FALSE)

rnamap <- data.frame(
    primary = c("Jack", "Jill", "Bob", "Barbara"),
    colname = c("samparray1", "samparray2", "samparray3",
                               "samparray4"),
    stringsAsFactors = FALSE)

## Combine as a named list and convert to a DataFrame
mylist <- list(Affy = exprmap, Methyl450k = methylmap,
               CNVgistic = rangemap, RNASeqGene = rnamap)

## Create a sampleMap
mySampleMap <- listToMap(mylist)
## Create an example phenotype data
pDat <- data.frame(sex = c("M", "F", "M", "F"),
                       age = 38:41,
                       row.names = c("Jack", "Jill", "Bob", "Barbara"))

## Create a MultiAssayExperiment instance
myMultiAssayExperiment <- MultiAssayExperiment(experiments = ExpList,
                                               pData = pDat,
                                               sampleMap = mySampleMap)
