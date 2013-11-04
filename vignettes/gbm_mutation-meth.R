
library(Biobase)
library(impute)
if(!require(copynumbR)){
    library(devtools)
    install_github("copynumbR", user="lima1")
    library(copynumbR)
}

##library(WGCNA)
## Install eigenR2 from http://www.genomine.org/eigenr2/
## library(eigenR2)
load("../inst/extdata/tcga_gbm/eset.sm.rda")
load("../inst/extdata/tcga_gbm/eset.meth27.rda")

##1 is silent, 2 is non-silent somatic mutations.  Set 1 to 0 (same as no mutation):
exprs(eset.sm)[exprs(eset.sm) == 1] <- 0

##discard genes with fewer than 5 somatic mutations - 76 left:
eset.sm <- eset.sm[, sampleNames(eset.sm) %in% sampleNames(eset.meth27)]
eset.sm <- eset.sm[esApply(eset.sm, 1, function(x) sum(x > 0) >= 5), ]

##Get rid of methylation features with fewer than 20 non-NA values (very bimodal):
eset.meth27 <- eset.meth27[esApply(eset.meth27, 1, function(x) sum(!is.na(x)) > 20), ]

samples.intersect <- intersect(sampleNames(eset.sm), sampleNames(eset.meth27))
features.intersect <- intersect(featureNames(eset.sm), featureNames(eset.meth27))
length(samples.intersect)  #116
length(features.intersect)  #39
eset.sm <- eset.sm[features.intersect, samples.intersect]
eset.meth27 <- eset.meth27[features.intersect, samples.intersect]

sm.meth.tt <- t(sapply(1:nrow(eset.sm), function(i){
    tmp <- t.test(exprs(eset.meth27[i, ]), exprs(eset.sm[i, ]))
    output <- c(tmp$p.value, diff(tmp$estimate))
    names(output) <- c("p.value", "diffmean")
    output
}))
rownames(sm.meth.tt) <- featureNames(eset.sm)

##A bunch are positively correlated (p<0.05 Bonferroni)
rownames(sm.meth.tt)[sm.meth.tt[, "p.value"] < 0.05/nrow(sm.meth.tt) & sm.meth.tt[, "diffmean"] < 0]

##But only EGFR is negatively correlated:
rownames(sm.meth.tt)[sm.meth.tt[, "p.value"] < 0.05/nrow(sm.meth.tt) & sm.meth.tt[, "diffmean"] > 0]

sm.meth.egfr <- rbind(exprs(eset.meth27)["EGFR", samples.intersect], exprs(eset.sm)["EGFR", ])
rownames(sm.meth.egfr) <- c("methylation", "mutation")
sm.meth.egfr.discrete <- sm.meth.egfr
sm.meth.egfr.discrete["methylation", ] <- ifelse(sm.meth.egfr.discrete[1, ] > 0.03, 2, 0)

pdf("gbm_EGFR_sm_meth_boxplot.pdf", width=3, height=3)
par(mar=c(3, 4, 0.5, 0.5))
boxplot(sm.meth.egfr["methylation", ] ~ sm.meth.egfr["mutation", ],
        ylim=c(0, 0.06), bw=0.1, names=c("WT/silent", "Mutated"), col="lightgrey",
        ylab="EGFR promoter Methylation")
dev.off()

pdf("gbm_EGFR_sm_meth_heatmap.pdf", width=6, height=6)
library(gplots)
heatmap.2(sm.meth.egfr.discrete, scale="none", col=c("lightgrey", "black"), trace="none",
          labCol=rep("", ncol(sm.meth.egfr.discrete)), key=FALSE, margins=c(1, 7), cexRow=1,
          rowsep=1:nrow(sm.meth.egfr.discrete), colsep=1:ncol(sm.meth.egfr.discrete), dendrogram="none")
par(xpd=NA)
legend(x=0.5, y=0.9, xjust=0.5, yjust=0.5, legend=c("non-synonymous mutation / hypermethylation", "silent or wild-type / normal methylation"), bty='n', pch=15, col=c("black", "lightgrey"))
dev.off()

## -------------------------------------------------
## EGFRvIII exon analysis
## -------------------------------------------------
load("../inst/extdata/tcga_gbm/eset.expr.agilent.all.rda")
load("../inst/extdata/tcga_gbm/eset.CN.gene.all.rda")
intersect.samples <- intersect(sampleNames(eset.CN.gene.all), sampleNames(eset.expr.agilent.all))
eset.CN.gene.all <- eset.CN.gene.all[, intersect.samples]
eset.expr.agilent.all <- eset.expr.agilent.all[, intersect.samples]
x <- exprs(eset.sm)["EGFR",]
eset.CN.gene.all$somatic <- c("Wildtype", "Silent","Non-Silent")[x[match(sampleNames(eset.CN.gene.all), names(x))]+1]

## Status of EGFRvIII mutation was determined as previously described
## (Szerlip et al., 2012). We compared the mean aCGH copy number
## ratios of probes in the deleted region of EGFR (Agilent probe ids
## A_14_P102368, A_16_P17952602, A_16_P38034515, A_16_P01718265), with
## the mean ratios of the probes in the 3’ prime end of EGFR,
## downstream of exon 7 (A_14_P133869, A_16_P17952748, A_16_P01718353,
## A_16_P17952840, A_16_P38034765, A_14_P106592). Difference of > 1
## standard deviation was considered significant, ie, exons 3-7 are
## deleted in these samples.

e3 <- read.csv("../inst/extdata/tcga_gbm/sm_exonsv2-v7.csv", header=TRUE, row.names=1, comment.char="#")
e3 <- na.omit(rownames(e3)[e3$EGFRvIII])
e33 <- intersect(e3, sampleNames(eset.CN.gene.all))
eset.CN.gene.all$somatic[match(e33, sampleNames(eset.CN.gene.all))] <- "EGFRvIII"

eset.CN.gene.all$somatic[is.na(eset.CN.gene.all$somatic)] <- "N/A"

eset.CN.gene.all$somatic[which(eset.CN.gene.all$somatic=="Wildtype") ] <- "Wildtype/Silent"
eset.CN.gene.all$somatic[which(eset.CN.gene.all$somatic=="Silent") ] <- "Wildtype/Silent"

obj <- copynumbR.boxplot(eset.CN.gene.all["EGFR",], eset.expr.agilent.all["EGFR",], ylab="Agilent mRNA Expression")

pdf("gbm_EGFRVIII.pdf", width=6, height=6)
library(ggplot2)
plot(obj$plot+geom_jitter(aes(color=eset.CN.gene.all[, id]$somatic, shape=eset.CN.gene.all[, id]$somatic), size=2)
     + scale_color_discrete(name="Mutation") + scale_shape_discrete(name="Mutation")
     + theme_classic2(16)+theme(axis.text.x= element_text(angle = 45,  hjust=1)))
dev.off()


## second question: which mutations affect global methylation?  This
## didn't really work out - skip for now.  Should have found IDH1,
## although only 5 samples with this mutation here.
##
## Next time try using eigenr2 (leek's website)
load("../inst/extdata/tcga_gbm/eset.meth27.rda")
meth.imputed <- impute.knn(exprs(eset.meth27))$data
meth.pc <- prcomp(t(meth.imputed))
scores.pc <- t(meth.pc$x)[, samples.intersect]


sm.methpc.tt <- lapply(1:10, function(j){  ##1:10 for PC1-10
    pc1vsmut <- t(sapply(1:nrow(eset.sm), function(i){
        tmp <- t.test(scores.pc[j, ], exprs(eset.sm[i, ]))
        output <- c(tmp$p.value, diff(tmp$estimate))
        names(output) <- c("p.value", "diffmean")
        output
    }))
    rownames(pc1vsmut) <- featureNames(eset.sm)
    pc1vsmut
})

mut.meth.shifts <- lapply(sm.methpc.tt, function(x){
    data.frame(x)[x[, "p.value"] < 0.05 / (10*nrow(x)), ]
})
names(mut.meth.shifts) <- paste("PC", 1:length(mut.meth.shifts), sep="")
mut.meth.shifts <- mut.meth.shifts[sapply(mut.meth.shifts, function(x) nrow(x) > 0)]
mut.meth.shifts  ##EGFR and TTN only - PC6, 7, 9.


hist(sm.meth.cor, breaks="FD")
interesting.genes <- c("EGFR", "PDGFRA", "NF1", "CDKN2A")
abline(v=sm.meth.cor[interesting.genes], col="red")

head(sort(sm.meth.cor), 15)
sm.meth.cor[interesting.genes]
