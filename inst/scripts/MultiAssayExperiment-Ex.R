## Run the example Elist
example("Elist")

## Load example GRangesList object
example("RangedRaggedAssay")

## Add the RangedRaggedAssay to the list
ExpList <- c(ExpList, myRRA)
names(ExpList)[3] <- "CNVgistic"

## Run the sample map example
example("sampleMap")

## Create an example phenotype data
myPheno <- data.frame(sex = c("M", "F", "M", "F"),
                       age = 38:41,
                       row.names = c("Jack", "Jill", "Bob", "Barbara"))

## Auto-creation of map
myMultiAssayExperiment <- MultiAssayExperiment(Elist = ExpList,
                                               pData = myPheno,
                                               sampleMap = mySampleMap)
