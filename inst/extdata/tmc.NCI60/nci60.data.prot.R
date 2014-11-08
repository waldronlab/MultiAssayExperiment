
# date: 05/05/2014
# aim: 
# - select proteins with log(Intensity,10) >= qnorm(0.01) in >= 5 cell lines
# - plot the average proteins log(Intensity,10) in the pnael 
#   by the frequency of non-N/A values 
# return: 
# - NCI60.prot.expr.pdf
# - NCI60.prot.txt
# - NCI60.prot.max5utr.txt
# - NCI60.prot.max3utr.txt


library("GenomicRanges")
library("gplots")
library("nortest")


load("summarizedExperiment.rda")
load("data.rda")
load("annotation.rda")

# proteomics data
exprm = as.matrix(assays(sset)[[2]])
# remove matrix rows for IDs artificially added 
# to costruct the summarizedExperiment object
ind = apply(exprm, 1, function(x) all(is.na(x)))
exprm = exprm[!ind,]

# transcriptomics data
rna.exprm = as.matrix(assays(sset)[[1]])

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


pdf("NCI60.prot.expr.pdf")

barplot(names.arg=names(table(nonNA)),height=as.numeric(table(nonNA)),
 xlab="Frequency of non N/A values per protein",cex.names=.5,las=2,
  ylab="N. proteins",col=c(rep(rgb(1,0,0),times=5),rep(rgb(1,1,1),times=54))) 

# remove IDs with fewer than 5 non-N/A values overall cell lines
exprm = exprm[nonNA>=5,]

# plot the relationship between frequency of non N/A values per protein
# and average protein log(Intensity,10)
nonNAcounts=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))
nonNAcountsTable=table(nonNAcounts)
out=list()
for(i in 1 : length(table(nonNAcounts))){
 tmp=names(nonNAcounts)[as.numeric(nonNAcounts)==
  names(nonNAcountsTable)[i]]
 tmp.exprm=exprm[match(tmp,rownames(exprm)),]
 out[[paste("nonNA_",names(nonNAcountsTable)[i],sep="")]]=
  as.numeric(apply(tmp.exprm,1,function(x) mean(x[is.na(x)==FALSE])))
 show(i)
}
boxplot(out,names=gsub("nonNA_","",names(out)),
 ylab="Mean protein log(Intensity,10)",cex=.75,pch=19,las=2,cex.lab=1,
  cex.axis=.55,main="No filter on protein intensity",
   xlab="Frequency of non-N/A values per protein")

# plot the relationship between frequency of non N/A values per protein
# and average transcript log(Intensity,10)
nonNAcounts=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))
nonNAcountsTable=table(nonNAcounts)
out=list()
for(i in 1 : length(table(nonNAcounts))){
 tmp=names(nonNAcounts)[as.numeric(nonNAcounts)==
  names(nonNAcountsTable)[i]]
 tmp.exprm=rna.exprm[match(tmp,rownames(rna.exprm)),]
 out[[paste("nonNA_",names(nonNAcountsTable)[i],sep="")]]=
  as.numeric(apply(tmp.exprm,1,function(x) mean(x[is.na(x)==FALSE])))
 show(i)
}
boxplot(out,names=gsub("nonNA_","",names(out)),
 ylab="Mean transcript log(Intensity,10)",cex=.75,pch=19,las=2,cex.lab=1,
  cex.axis=.55,main="No filter on protein intensity", 
   xlab="Frequency of non-N/A values per protein")


# plot log(Intensity,10) by cell line before filtering on intesity
tmp=matrix(nrow=dim(exprm)[1],ncol=dim(exprm)[2],byrow=TRUE)
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(is.na(exprm[i,j])==TRUE){tmp[i,j]=0}
  if(is.na(exprm[i,j])==FALSE){tmp[i,j]=exprm[i,j]}
 }
}
heatmap.2(tmp,trace="none",col=colorRampPalette(c(rgb(1,1,1),
 rgb(0,1,0),rgb(0,0,1))),labRow=rep("",dim(exprm)[1]),
  main="No filter on protein intensity",labCol=rep("",dim(exprm)[2]),
   dendro="row",xlab="Cell line",ylab="Protein")


# test for normality of log(Intensity,10) distribution
cvm.test(as.numeric(exprm))
summary(as.numeric(exprm))
SD=sd(as.numeric(exprm)[is.na(as.numeric(exprm)) == FALSE])
AVG=mean(as.numeric(exprm)[is.na(as.numeric(exprm)) == FALSE])
# plot log(Intensity,10) distribution
hist(as.numeric(exprm),breaks=25,las=2,cex.lab=.75,prob=TRUE,
 xlab="Protein log(Intensity,10)",main="")
abline(v=qnorm(0.01,sd=SD,mean=AVG),lty=2,col=rgb(1,0,0))
# remove IDs with log(intensity,10) values <=qnorm(0.01) in >=5 cell lines
NumLineLowE=apply(exprm,1,function(x) length(x[x<=qnorm(0.01,sd=SD,mean=AVG)]))
exprm=exprm[NumLineLowE<5,]
# plot log(Intensity,10) by cell line after filter
tmp=matrix(nrow=dim(exprm)[1],ncol=dim(exprm)[2],byrow=TRUE)
for(i in 1 : dim(exprm)[1]){
 for(j in 1 : dim(exprm)[2]){
  if(is.na(exprm[i,j])==TRUE){tmp[i,j]=0}
  if(is.na(exprm[i,j])==FALSE){tmp[i,j]=exprm[i,j]}
 }
}
heatmap.2(tmp,trace="none",na.rm=TRUE,col=colorRampPalette(c(rgb(1,1,1),
 rgb(0,1,0),rgb(0,0,1))),labRow=rep("",dim(exprm)[1]),
  main="Filter on protein intensity",labCol=rep("",dim(exprm)[2]),
   dendro="row",xlab="Cell line",ylab="Protein")


# plot the relationship between frequency of non N/A values per protein 
# and average proteins log(Intensity,10) 
nonNAcounts=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))
nonNAcountsTable=table(nonNAcounts)
out=list()
for(i in 1 : length(table(nonNAcounts))){
 tmp=names(nonNAcounts)[as.numeric(nonNAcounts)==
  names(nonNAcountsTable)[i]]
 tmp.exprm=exprm[match(tmp,rownames(exprm)),]
 out[[paste("nonNA_",names(nonNAcountsTable)[i],sep="")]]=
  as.numeric(apply(tmp.exprm,1,function(x) mean(x[is.na(x)==FALSE])))
 show(i) 
}
boxplot(out,names=gsub("nonNA_","",names(out)),
 ylab="Mean protein log(Intensity,10)",cex=.75,pch=19,las=2,cex.lab=1,
  cex.axis=1,main="Filter on protein intensity",
   xlab="Frequency of non-N/A values per protein")

# plot the relationship between frequency of non N/A values per protein
# and average proteins log(Intensity,10)
nonNAcounts=apply(exprm,1,function(x) length(x[is.na(x)==FALSE]))
nonNAcountsTable=table(nonNAcounts)
out=list()
for(i in 1 : length(table(nonNAcounts))){
 tmp=names(nonNAcounts)[as.numeric(nonNAcounts)==
  names(nonNAcountsTable)[i]]
 tmp.exprm=rna.exprm[match(tmp,rownames(rna.exprm)),]
 out[[paste("nonNA_",names(nonNAcountsTable)[i],sep="")]]=
  as.numeric(apply(tmp.exprm,1,function(x) mean(x[is.na(x)==FALSE])))
 show(i)
}
boxplot(out,names=gsub("nonNA_","",names(out)),
 ylab="Mean transcript log(Intensity,10)",cex=.75,pch=19,las=2,cex.lab=1,
  cex.axis=1,main="Filter on protein intensity",
   xlab="Frequency of non-N/A values per protein")

dev.off()


# print selected proteins
prots=unique(annotation$HGNC.symbol[match(rownames(exprm),annotation$ID)])
write.table(prots,"NCI60.prot.txt",row.names=FALSE,
 col.names=FALSE,quote=FALSE)

# print UCSC ids of maximal 5' UTR corresponding to selected proteins
prot.ucsc.max5utr=unique(annotation$UCSC.max5utr[match(rownames(exprm),
 annotation$ID)])
prot.ucsc.max5utr=prot.ucsc.max5utr[is.na(prot.ucsc.max5utr)==FALSE]
write.table(prot.ucsc.max5utr,"NCI60.prot.max5utr.txt",
 row.names=FALSE,col.names=FALSE,quote=FALSE)

# print UCSC ids of maximal 3' UTR corresponding to selected proteins
prot.ucsc.max3utr=unique(annotation$UCSC.max3utr[match(rownames(exprm),
 annotation$ID)])
prot.ucsc.max3utr=prot.ucsc.max3utr[is.na(prot.ucsc.max3utr)==FALSE]
write.table(prot.ucsc.max3utr,"NCI60.prot.max3utr.txt",
 row.names=FALSE,col.names=FALSE,quote=FALSE)



