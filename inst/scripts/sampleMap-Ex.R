## Create sample maps for each experiment
exprmap <- data.frame(
    master = c("Jack", "Jill", "Barbara", "Bob"),
    assay = c("array1", "array2", "array3", "array4"), 
    stringsAsFactors = FALSE)
methylmap <- data.frame(
    master = c("Jack", "Jack", "Jill", "Barbara", "Bob"),
    assay = c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5"),
    stringsAsFactors = FALSE)
rangemap <- data.frame(master = c("Jack", "Jill", "Jill"), 
    assay = c("snparray1", "snparray2", "snparray3"),
    stringsAsFactors = FALSE)

## Combine as a named list and convert to a DataFrame
mylist <- list(exprmap, methylmap, rangemap)
names(mylist) <- c("Affy", "Methyl450k", "CNVgistic")

## Create a sampleMap
mySampleMap <- listToMap(mylist)
