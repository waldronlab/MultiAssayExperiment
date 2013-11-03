library(Biobase)
library(sva)

load("../inst/extdata/tcga_ov/ovc_eset_copynumber_gistic_2013_01_16.rda")
load("../inst/extdata/tcga_ov/TCGA_eset.rda")
##load("../inst/extdata/tcga_ov/TCGA.RNASeqV2_eset.rda")
load("../inst/extdata/tcga_ov/TCGA.agilent_eset.rda")
source("../inst/corFinder.R")

all.samples <- intersect(sampleNames(TCGA_eset), intersect(sampleNames(TCGA.agilent_eset), sampleNames(eset.cn)))
all.features <- intersect(featureNames(TCGA_eset), intersect(featureNames(TCGA.agilent_eset), featureNames(eset.cn)))
TCGA_eset <- TCGA_eset[all.features, all.samples]
TCGA.agilent_eset <- TCGA.agilent_eset[all.features, all.samples]
eset.cn <- eset.cn[all.features, all.samples]

cn.cor <- corFinder(list(microarray=TCGA_eset, cn=eset.cn), use.ComBat=FALSE)
##rnaseq.cor <- corFinder(list(rnaseq=TCGA.RNASeqV2_eset, microarray=TCGA_eset), use.ComBat=TRUE)
agilent.cor <- corFinder(list(microarray=TCGA_eset, agilent=TCGA.agilent_eset), use.ComBat=TRUE)


pdf("copy_number_QC.pdf", width=3, height=3)
par(mar=c(4, 4, 0.1, 0.1))
hist(diag(cn.cor), breaks="FD", main="", xlab="correlation")
abline(v=0.066, col="red", lw=2)
arrows(x0=0.06, y0=9, x1=0.002, y1=9, col="red", lw=2, length=0.1)
text(x=0.025, y=c(25, 18), labels=c("flagged", "samples"), col="red")
dev.off()

hist(diag(agilent.cor), breaks="FD", xlim=c(0.59, 1))
abline(v=0.87, col="red", lw=2)

sum(diag(cn.cor) < 0.066)
sum(diag(cn.cor) < 0.066 & diag(agilent.cor) < 0.87)

qqplot(diag(cn.cor), diag(agilent.cor), ylim=c(0.59, 1), pch=".")
abline(v=0.066, col="red", lw=2)

expronly <- corFinder(list(TCGA_eset, TCGA_eset), use.ComBat=FALSE)
hist(expronly[upper.tri(expronly)], xlim=c(0.8, 1), breaks="FD")
summary(as.numeric(expronly) > 0.98)
