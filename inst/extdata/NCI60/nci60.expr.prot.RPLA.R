
# Last updated on 10/12/2014


require(affy)
require(biocMultiAssay)
require(data.table)


##Recipe for preparing NCI60 text file

file.name <- "nci60_Protein__Lysate_Array_log2.txt.zip"
outfile.name <- "nci60_prot_RPLA_log2.csv"

url <- paste("http://discover.nci.nih.gov/cellminerdata/normalizedArchives/", file.name, sep="")
zipfile.name <- tempfile()
temp.dir <- tempdir()

download.file(url, destfile=zipfile.name)
unzip(zipfile.name, exdir=temp.dir)

##input filename:
textfile.name <- file.path(temp.dir, sub("\\.zip", "", file.name))

##column 1 is all identical entries, 3 and 4 are gene symbols and entrez IDs, so drop them.
affydat <- fread(textfile.name)
rpladat <- fread(textfile.name, drop=c(1, 3, 4), na.strings="-")
##one column of all NAs:
keep.cols <- sapply(rpladat, function(x) sum(!is.na(x))) > 0
rpladat <- rpladat[, keep.cols, with=FALSE]
setkey(rpladat, V2)
##average rows with identical probeset IDs:
rpladat <- rpladat[, lapply(.SD, mean), by=V2]

##now set column names, which fread misses:
cnames <- readLines(textfile.name, n=1)
cnames <- c("V2", strsplit(cnames, split="\t")[[1]][-1:-4])
cnames <- cnames[keep.cols]
cnames[cnames=="V2"] <- "id"
setnames(rpladat, colnames(rpladat), cnames)

##Create SummarizedExperiment
library(affy)
eset <- ExpressionSet(assayData=as.matrix(rpladat))
featureNames(eset) <- rpladat$id
annotation(eset) <- "hgu133plus2"
#library(biocMultiAssay)
#nci60.expr.se <- exs2se(eset)

##Write matrix and SE to file
write.csv(rpladat, row.names=FALSE, file=outfile.name)
#save(nci60.expr.se, file="nci60.expr.se.rda")

