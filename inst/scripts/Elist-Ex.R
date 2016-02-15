## Create an empty Elist instance 
Elist()

## Create array matrix and AnnotatedDataFrame to create an ExpressionSet class
arraydat <- matrix(seq(101, 108), ncol=4,
                   dimnames = list(
                     c("ENST00000294241", "ENST00000355076"), 
                     c("array1", "array2", "array3", "array4")
                   ))
arraypdat <- as(data.frame(
  slope53 = rnorm(4), 
  row.names = c("array1", "array2", "array3", "array4")), 
  "AnnotatedDataFrame")
exprdat <- Biobase::ExpressionSet(assayData=arraydat, phenoData=arraypdat)

## Create a sample methylation dataset
methyldat <- matrix(1:10,
                    ncol = 5, 
                    dimnames = list(
                      c("ENST00000355076", "ENST00000383706"),
                      c("methyl1", "methyl2", "methyl3", "methyl4", "methyl5")))

## Combine to a named list and call the Elist constructor function
ExpList<- list(arraydat, methyldat)
names(ExpList) <- c("Affy", "Methyl450k")
myElist <- Elist(ExpList)
