
# date: 14/05/2014
# aim: 
# - plot histogram of RMSPE ratios, M0 RMSPE, M1 RMSPE
# - plot histogram of correlation between preditions and 
#   observations achieved by M1 in testing subset
# - plot histogram of correlation between preditions and
#   observations achieved by M0 in testing subset
# - plot histogram of correlation ratios
#
# return: 
# - extra.sample.predict.out.pdf


library("GenomicRanges")
library("gplots")


load("../../results/summarizedExperiment.rda")
load("../../results/data.rda")
load("../../results/annotation.rda")


# proteomics data
exprm = as.matrix(assays(sset)[[2]])
# remove matrix rows for IDs artificially added
# to costruct the summarizedExperiment object
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]

# replace 0 values with N/A according to the specification by Prof. Kuster
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(exprm[i,j]==0){ exprm[i,j]=NA }
 }
}

# remove IDs whose protein levels are N/A overall cell lines
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]

# count non-N/A values per gene overall cell lines
nonNA=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))


# summary of extra-samples prediction error performances 
cvf = read.table("../../results/model/repeated.5fold.cv.conf.int/fit.summary.repeated.5fold.cv.knn.txt",header=T,sep="\t",as.is=TRUE)


pdf("extra.sample.predict.out.pdf")
par(mfrow=c(2,3))
# plot histogram of average RMSPE ratios, M0 RMSPE, M1 RMSPE
tmp = cbind(cvf$HGNC.symbol,cvf$rmspe.ratio.avg)
tmp = unique(tmp)
rmspe.ratio.avg=as.numeric(tmp[,2])
hist(rmspe.ratio.avg,breaks=25,prob=T,main="",xlab="Average M0/M1 RMSPE ratio")

tmp = cbind(cvf$HGNC.symbol,cvf$m0.rmspe.avg)
tmp = unique(tmp)
m0.rmspe.avg=as.numeric(tmp[,2])
hist(m0.rmspe.avg,breaks=25,prob=T,main="",xlab="Average M0 RMSPE")

tmp = cbind(cvf$HGNC.symbol,cvf$m1.rmspe.avg)
tmp = unique(tmp)
m1.rmspe.avg=as.numeric(tmp[,2])
hist(m1.rmspe.avg,breaks=25,prob=T,main="",xlab="Average M1 RMSPE")

# histogram of average ratios of correlations between pred and obs 
tmp = cbind(cvf$HGNC.symbol,cvf$cor.ratio.avg)
tmp = unique(tmp)
cor.ratio.avg=as.numeric(tmp[,2])
hist(cor.ratio.avg,breaks=100,prob=T,main="",
 xlab="Average M1/M0 cor(predict,obs) ratio")

tmp = cbind(cvf$HGNC.symbol,cvf$m0.cor)
tmp = unique(tmp)
m0.cor=as.numeric(tmp[,2])
hist(m0.cor,breaks=25,prob=T,main="",xlab="M0 cor(predict,obs)")

tmp = cbind(cvf$HGNC.symbol,cvf$X.m1.cor)
tmp = unique(tmp)
m1.cor=as.numeric(tmp[,2])
hist(m1.cor,breaks=25,prob=T,main="",xlab="M1 cor(predict,obs)")

dev.off()










