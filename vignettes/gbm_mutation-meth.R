library(Biobase)
library(impute)
library(WGCNA)



##load("../inst/extdata/tcga_gbm/eset.CN.region.all.rda")
load("../inst/extdata/tcga_gbm/eset.sm.rda")
load("../inst/extdata/tcga_gbm/eset.meth27.rda")

##1 is silent, 2 is non-silent somatic mutations.  Set 1 to 0 (same as no mutation):
exprs(eset.sm)[exprs(eset.sm) == 1] <- 0

##discard genes with fewer than 10 somatic mutations - only 39 left:
eset.sm <- eset.sm[esApply(eset.sm, 1, function(x) sum(x) >= 10), ]

##Get rid of methylation features with fewer than 20 non-NA values (very bimodal):
eset.meth27 <- eset.meth27[esApply(eset.meth27, 1, function(x) sum(!is.na(x)) > 20), ]
##exprs(eset.meth27) <- impute.knn(exprs(eset.meth27))$data

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

##A bunch are positively correlated:
rownames(sm.meth.tt)[sm.meth.tt[, "p.value"] < 1e-3 & sm.meth.tt[, "diffmean"] < 0]
## [1] "C3"      "CDH9"    "CNTNAP2" "COL6A3"  "DRD5"    "FBN3"    "GABRA6"
## [8] "HSPG2"   "MUC17"   "MUC5B"   "OBSCN"   "PIK3R1"  "PKHD1"   "PPP1R3A"
##[15] "RYR2"    "SCN9A"   "TRPV6"   "USH2A"   "VWF"

##But only EGFR is negatively correlated:
rownames(sm.meth.tt)[sm.meth.tt[, "p.value"] < 1e-3 & sm.meth.tt[, "diffmean"] > 0]

## second question:  which mutations affect global methylation?
## Use eigenr2 (leek's website)
meth.pc <- prcomp(exprs(eset.meth27))




for (i in which(sm.meth.tt[, "p.value"] < 1e-3 & sm.meth.tt[, "diffmean"] > 0)){
    tmp <- data.frame(methylation=exprs(eset.meth27)[i, ], sm=factor(exprs(eset.sm)[i, ]))
    boxplot(methylation ~ sm, data=tmp, main=featureNames(eset.sm)[i])
}

hist(sm.meth.cor, breaks="FD")
interesting.genes <- c("EGFR", "PDGFRA", "NF1", "CDKN2A")
abline(v=sm.meth.cor[interesting.genes], col="red")

head(sort(sm.meth.cor), 15)
sm.meth.cor[interesting.genes]

tmp <- data.frame(sm=factor(exprs(eset.sm["TRAPPC9", ])), meth=exprs(eset.meth27["TRAPPC9", ]))
