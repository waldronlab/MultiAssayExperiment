## Create an empty ExperimentList instance
ExperimentList()

## Create array matrix and AnnotatedDataFrame to create an ExpressionSet class
arraydat <- matrix(data = seq(101, length.out = 20), ncol = 4,
    dimnames = list(
        c("ENST00000294241", "ENST00000355076",
        "ENST00000383706","ENST00000234812", "ENST00000383323"),
        c("array1", "array2", "array3", "array4")
    ))

colDat <- data.frame(slope53 = rnorm(4),
    row.names = c("array1", "array2", "array3", "array4"))

## SummarizedExperiment constructor
exprdat <- SummarizedExperiment::SummarizedExperiment(arraydat,
    colData = colDat)

## Create a sample methylation dataset
methyldat <- matrix(data = seq(1, length.out = 25), ncol = 5,
    dimnames = list(
        c("ENST00000355076", "ENST00000383706",
          "ENST00000383323", "ENST00000234812", "ENST00000294241"),
        c("methyl1", "methyl2", "methyl3",
          "methyl4", "methyl5")
    ))

## Create a sample RNASeqGene dataset
rnadat <- matrix(
    data = sample(c(46851, 5, 19, 13, 2197, 507,
        84318, 126, 17, 21, 23979, 614), size = 20, replace = TRUE),
    ncol = 4,
    dimnames = list(
        c("XIST", "RPS4Y1", "KDM5D", "ENST00000383323", "ENST00000234812"),
        c("samparray1", "samparray2", "samparray3", "samparray4")
    ))

## Create a mock RangedSummarizedExperiment from a data.frame
rangedat <- data.frame(chr="chr2", start = 11:15, end = 12:16,
    strand = c("+", "-", "+", "*", "."),
    samp0 = c(0,0,1,1,1), samp1 = c(1,0,1,0,1), samp2 = c(0,1,0,1,0),
    row.names = c(paste0("ENST", "00000", 135411:135414), "ENST00000383323"))

rangeSE <- SummarizedExperiment::makeSummarizedExperimentFromDataFrame(rangedat)

## Combine to a named list and call the ExperimentList constructor function
assayList <- list(Affy = exprdat, Methyl450k = methyldat, RNASeqGene = rnadat,
                GISTIC = rangeSE)

## Use the ExperimentList constructor
ExpList <- ExperimentList(assayList)
