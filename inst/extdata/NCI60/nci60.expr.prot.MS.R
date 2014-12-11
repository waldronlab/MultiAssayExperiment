
##Retrieve NCI-60 MS-based protein data summarized by gene in CSV file 
##from http://129.187.44.58:7070/NCI60/main/download
##Three data types were available and are returned separately in CSV files:
##proteomes, deep proteomes and kinomes


require(affy)
require(biocMultiAssay)
require(data.table)


##input NCI-60 MS-based protein data
load("nci60_prot_MS_log10.rda")

summarize <- function(x,y,z){
 cnames <- y
 dat <- x
 outfile.name <- z
##select protein data quantified by iBAQ approach
cnames <- gsub(" ",".",cnames)
keep.cols <- c("Fasta.headers",cnames[grepl("iBAQ",cnames)==T])
dat <- dat[, match(keep.cols,cnames)]

ids <- unlist(lapply(dat[,1],
 function(x)(strsplit(strsplit(x,"Gene_Symbol=")[[1]][2]," ")[[1]][1])))
dat[,1] <- ids

##remove entries without valid HGNC symbols 
dat <- dat[is.na(dat[,1])==F,]
dat <- dat[is.na(match(dat[,1],"-"))==T,]
##null values correspond to undetected proteins 
dat[dat==0] <- NA

dat <- data.table(dat)
setkey(dat, "Fasta.headers")
##average rows with identical HGNC symbols:
dat <- dat[, lapply(.SD, mean), by="Fasta.headers"]

setnames(dat, colnames(dat), gsub("Fasta.headers","id",keep.cols))

##Create ExpressionSet
eset <- ExpressionSet(assayData=as.matrix(dat))
featureNames(eset) <- dat$id

##Write matrix 
write.csv(dat, row.names=FALSE, file=outfile.name)
}

summarize(dat,cnames,"nci60_prot_MS_log10.csv")
summarize(deepdat,deepcnames,"nci60_deep_prot_MS_log10.csv")
summarize(kindat,cnames,"nci60_kin_prot_MS_log10.csv")


