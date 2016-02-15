## Run the example Elist
example("Elist-class")

## Run the sample map example
example("listToMap")

## Create an example phenotype data
myPheno <- data.frame(sex = c("M", "F", "M", "F"),
                       age = 38:41,
                       row.names = c("Jack", "Jill", "Bob", "Barbara"))

## Auto-creation of map
myMultiAssayExperiment <- MultiAssayExperiment(Elist = myElist,
                                               pData = myPheno,
                                               sampleMap = myDF)
