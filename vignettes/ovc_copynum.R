library(Biobase)
library(sva)
library(impute)

load("../inst/extdata/tcga_ov/ovc_eset_copynumber_gistic_2013_01_16.rda")
load("../inst/extdata/tcga_ov/TCGA_eset.rda")
##load("../inst/extdata/tcga_ov/TCGA.RNASeqV2_eset.rda")
load("../inst/extdata/tcga_ov/TCGA.agilent_eset.rda")
source("../inst/corFinder.R")

cn.cor <- corFinder(list(microarray=eset.expression, cn=eset.cn), use.ComBat=FALSE)
##rnaseq.cor <- corFinder(list(rnaseq=TCGA.RNASeqV2_eset, microarray=TCGA_eset), use.ComBat=TRUE)
agilent.cor <- corFinder(list(microarray=TCGA_eset, agilent=TCGA.agilent_eset), use.ComBat=TRUE)
colnames(cn.cor) <- sub(".+:", "", colnames(cn.cor))
rownames(cn.cor) <- sub(".+:", "", rownames(cn.cor))
colnames(agilent.cor) <- sub(".+:", "", colnames(agilent.cor))
rownames(agilent.cor) <- sub(".+:", "", rownames(agilent.cor))

pdf("copy_number_QC.pdf", width=3, height=3)
par(mar=c(4, 4, 0.1, 0.1))
hist(diag(cn.cor), breaks="FD", main="", xlab="correlation")
abline(v=0.066, col="red", lw=2)
arrows(x0=0.06, y0=9, x1=0.002, y1=9, col="red", lw=2, length=0.1)
text(x=0.025, y=c(25, 18), labels=c("flagged", "samples"), col="red")
dev.off()

## Correlation of Affy - Agilent technical replicates
hist(diag(agilent.cor), breaks="FD", xlim=c(0.59, 1))
abline(v=0.87, col="red", lw=2)

## Counts of outliers for CN and for technical reps (qualitative cutoff)
samples.intersect <- intersect(colnames(cn.cor), colnames(agilent.cor))
sum(diag(cn.cor[samples.intersect, samples.intersect]) < 0.066)
sum(diag(cn.cor[samples.intersect, samples.intersect]) < 0.066 & diag(agilent.cor[samples.intersect, samples.intersect]) < 0.87)

## Methylation is not nearly so informative for QC in this case, but I
## left this analysis in place anyways.  It still identifies three
## outliers, two of which overlap with the CN outliers.  There is
## overall a strong negative correlation between expression and
## methylation (mean cor ~ -0.35)
load("../inst/extdata/tcga_ov/ovc_eset_methylation.rda")
eset.meth <- eset.meth[esApply(eset.meth, 1, function(x) sum(is.na(x))) < 300, ]
exprs(eset.meth) <- impute.knn(exprs(eset.meth))$data

meth.cor <- corFinder(list(microarray=TCGA_eset, meth=eset.meth), use.ComBat=FALSE)
colnames(meth.cor) <- sub(".+:", "", colnames(meth.cor))
rownames(meth.cor) <- sub(".+:", "", rownames(meth.cor))
hist(diag(meth.cor), breaks="FD")
sum(diag(meth.cor) > -0.2)

## Overlap between methylation outliers and copy number outliers.
samples.intersect <- intersect(colnames(meth.cor), colnames(cn.cor))
length(samples.intersect)
summary(diag(meth.cor[samples.intersect, samples.intersect]) > -0.2)
summary(diag(cn.cor[samples.intersect, samples.intersect]) < 0.066)
summary(diag(meth.cor[samples.intersect, samples.intersect]) > -0.2 & diag(cn.cor[samples.intersect, samples.intersect]) < 0.066)

