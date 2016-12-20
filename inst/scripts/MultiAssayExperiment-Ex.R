## Run the example ExperimentList
example("ExperimentList")

## Load example GRangesList object
example("RangedRaggedAssay")

## Add the RangedRaggedAssay to the list
ExpList <- c(ExpList, myRRA)
names(ExpList)[4] <- "CNVgistic"

## Run the sample map example
example("sampleMap")

## Create an example phenotype data
pDat <- data.frame(sex = c("M", "F", "M", "F"),
                       age = 38:41,
                       row.names = c("Jack", "Jill", "Bob", "Barbara"))

## Create a MultiAssayExperiment instance
myMultiAssayExperiment <- MultiAssayExperiment(experiments = ExpList,
                                               pData = pDat,
                                               sampleMap = mySampleMap)
