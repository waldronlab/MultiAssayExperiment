##Recipe for preparing NCI60 text file

file.name <- "nci60_RNA__Affy_HG_U133(A_B)_GCRMA.txt.zip"
outfile.name <- "affy_gcrma.csv"

url <- paste("http://discover.nci.nih.gov/cellminerdata/normalizedArchives/", file.name, sep="")
zipfile.name <- tempfile()
temp.dir <- tempdir()

download.file(url, destfile=zipfile.name)
unzip(zipfile.name, exdir=temp.dir)

##input filename:
textfile.name <- file.path(temp.dir, sub("\\.zip", "", file.name))

##column 1 is all identical entries, 3 and 4 are gene symbols and entrez IDs, so drop them.
library(data.table)
affydat <- fread(textfile.name, drop=c(1, 3, 4), na.strings="-")
##one column of all NAs:
keep.cols <- sapply(affydat, function(x) sum(!is.na(x))) > 0
affydat <- affydat[, keep.cols, with=FALSE]
setkey(affydat, V2)
##average rows with identical probeset IDs:
affydat <- affydat[, lapply(.SD, mean), by=V2]

##now set column names, which fread misses:
cnames <- readLines(textfile.name, n=1)
cnames <- c("V2", strsplit(cnames, split="\t")[[1]][-1:-4])
cnames <- cnames[keep.cols]
cnames[cnames=="V2"] <- "id"
setnames(affydat, colnames(affydat), cnames)

##Create SummarizedExperiment
library(affy)
eset <- ExpressionSet(assayData=as.matrix(affydat))
featureNames(eset) <- affydat$id
annotation(eset) <- "hgu133plus2"
library(biocMultiAssay)
nci60.expr.se <- exs2se(eset)

##Write matrix and SE to file
write.csv(affydat, row.names=FALSE, file=outfile.name)
save(nci60.expr.se, file="nci60.expr.se.rda")
