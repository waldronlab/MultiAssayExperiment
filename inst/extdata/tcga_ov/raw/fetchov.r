# Script for downloading Ovarian Cancer ("OV") TCGA data using RTCGAToolbox
if( !require(RTCGAToolbox) ){
  BiocInstaller::biocLite("LiNk-NY/RTCGAToolbox")
  library(RTCGAToolbox)
}
# Set download directory
mydir <- file.path("/scratch/")

# Read in data from RTCGAToolbox
source("tcgarundates.r")
load("rundates.rda")
recentrun <- gsub("-", "", dates$run[nrow(dates)])
recentgist <- gsub("-", "", dates$analyze[nrow(dates)])
ovres <- getFirehoseData("OV", runDate=recentrun, gistic2_Date=recentgist,
                          RNAseq_Gene=TRUE, Clinic=TRUE,
                          RNAseq2_Gene_Norm=TRUE, CNA_SNP=TRUE, CNV_SNP=TRUE, 
                          CNA_Seq=TRUE, CNA_CGH=TRUE,  Methylation=TRUE, 
                          Mutation=TRUE, mRNA_Array=TRUE, miRNA_Array=TRUE, RPPA=TRUE, destdir = mydir)

OVmirna <- getFirehoseData("OV", runDate=recentrun, miRNASeq_Gene=TRUE)

library(data.table)

grep("OV", dir(mydir), ignore.case=TRUE, value=TRUE)


Clin <- fread(file.path(mydir, paste0(recentrun, "-OV-Clinical.txt")), stringsAsFactors=FALSE)

for (i in 1:3) {
  Clin <- rbindlist(list(Clin, lapply(strsplit(colnames(Clin), "-"), "[", i)))
}

setnames(Clin, "Hybridization REF", "Information")

Clin[[1]][28:30] <- c("colsite", "tss", "partID")

# Code for renaming columns from TCGA barcode to n()
setnames(Clin, 2:length(Clin), as.character(seq(1:(length(Clin)-1))))

# number missing
sum(is.na(Clin))

# Read in OV All Thresholded by Genes File
AllGenes <- fread(file.path(mydir, paste0(recentgist, "-OV-all_thresholded.by_genes.txt")), stringsAsFactors=FALSE)

for (i in 1:3) {
  AllGenes <- rbindlist(list(AllGenes, lapply(strsplit(colnames(AllGenes), "-"), "[", i)))
}

AllGenes[[1]][(nrow(AllGenes)-2):nrow(AllGenes)] <- c("colsite", "tss", "partID")
AllGenes[[2]][(nrow(AllGenes)-2):nrow(AllGenes)] <- NA
AllGenes[[3]][(nrow(AllGenes)-2):nrow(AllGenes)] <- NA

# Code for renaming columns from TCGA barcode to n()
setnames(AllGenes, colnames(AllGenes)[grep("[0-9]", colnames(AllGenes))], as.character(seq(4:(length(AllGenes)))))

sum(is.na(AllGenes))

# Read Methylation File
Methyl <- fread(file.path(mydir, paste0(recentrun, "-OV-Methylation-1.txt")), stringsAsFactors=FALSE)

for (i in 1:3){
  Methyl <- rbindlist(list(Methyl, lapply(strsplit(colnames(Methyl), "-") , "[", i)))
}

Methyl[[1]][(nrow(Methyl)-2):nrow(Methyl)] <- c("colsite", "tss", "partID")

sum(is.na(Methyl))

# setnames(Methyl, colnames(Methyl)[grep("TCGA", colnames(Methyl))], as.character(seq(2:(length(Methyl)))))

###
# Duplicate Barcodes do not allow easy merging #
###
# library(reshape2)
# melt(Methyl, id.vars = "Hybridization REF")
# setnames(Methyl, colnames(Methyl)[-1], Methyl[1, 2:length(Methyl), with=FALSE])

# Data with fist duplicates
Methyl[, seq(1, 2445, 4), with=FALSE]

gg <- 1:2449  
gg <- gg[-c(seq(1,2445, 4))]

# Data with remaining duplicates
Methyl[, gg, with=FALSE]

# example
ext <- Methyl[1:10, 1:10, with=FALSE]
save(ext, file=file.path(mydir, "methylEX.rda"))

# Load OV RNAseqGene 

ovgist <- fread(file.path(mydir, paste0(recentrun, "-OV-RNAseqGene.txt")), stringsAsFactors=FALSE)


