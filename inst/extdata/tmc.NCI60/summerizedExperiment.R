
# date: 26/03/2014
# aim: define SummarizedExperiment-class object for NCI60 processed data. 
#  Metadata are RNA/protein gene expression values in log10 scale. 
# output: 
# - ./results/NCI60.summarizedExperiment.rda
# - ./results/NCI60.rna.expr.pdf
# - ./results/NCI60.prot.expr.pdf


library(GenomicRanges)
library(gplots)


load("data.rda")
load("annotation.rda")


# get gene identifiers in either RNA or protein assay
g.rna = annotation$ID[is.na(match(annotation$HGNC.symbol,rna[,3]))==FALSE]
g.prot = annotation$ID[is.na(match(annotation$IPI,prot$Accession))==FALSE]
g = union(g.rna,g.prot)
annotation = annotation[is.na(match(annotation$ID,g))==FALSE,]
# sort annotation by the order of genes in g
annotation = annotation[match(g, annotation$ID),]


# summarize probe-level intensity values (log10) to gene-level RNA intensity values
expr.rna = data.frame()
for(i in 1 : length(g)){
 # columns 9-67 correspond to NCI60 cell lines
 tmp = rna[is.na(match(rna[,3],annotation$HGNC.symbol[is.na(match(annotation$ID,g[i]))==FALSE]))==FALSE,9:67] 
 if(dim(tmp)[1]==0){ expr.rna = rbind(expr.rna, rep(NA,length(colnames(rna)[9:67]))) }
 if(dim(tmp)[1]!=0){ expr.rna = rbind(expr.rna,colMeans(tmp)) }
 show(i)
}
colnames(expr.rna) = colnames(rna)[9:67]
rownames(expr.rna) = g


# get (log10) protein values
expr.prot = data.frame() 
for(i in 1 : length(g)){
 tmp = prot[is.na(match(prot$Accession,annotation$IPI[is.na(match(annotation$ID,g[i]))==FALSE]))==FALSE,
  grepl("LFQ",colnames(prot))==TRUE]
 if(dim(tmp)[1]==0){ expr.prot = rbind(expr.prot, rep(NA,length(colnames(prot)[grepl("LFQ",colnames(prot))==TRUE]))) } 
 if(dim(tmp)[1]!=0){ expr.prot = rbind(expr.prot,colMeans(tmp)) }
 show(i)
}
colnames(expr.prot) = colnames(prot)[grepl("LFQ",colnames(prot))==TRUE]
colnames(expr.prot) = gsub("LFQ_","",colnames(expr.prot))
rownames(expr.prot) = g


# make summarizedExperiment object
rowData <- GRanges(seqnames = annotation$Chromosome.Name,ranges = IRanges(start=annotation$Gene.Start..bp.,
 end = annotation$Gene.End..bp., names = annotation$ID),strand = annotation$Strand)
colData <- DataFrame(row.names=colnames(expr.rna),ID=colnames(expr.rna))
sset <- SummarizedExperiment(assays=list(as.matrix(expr.rna,row.names=annotation$ID),
 as.matrix(expr.prot,row.names=annotation$ID)),rowData=rowData,colData=colData)

save(sset,file="summarizedExperiment.rda")



