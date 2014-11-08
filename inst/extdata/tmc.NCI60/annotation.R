
### date: 19/03/2014
### aim: 
### - retrieve IPI-HGNC relations directly from NCI60 files with RNA and protein processed data
### - retrieve genomic ccordinates for genes profiled at either level in the NCI60 panel
### output: 
### - annotation.rda


# load Ensembl BioMart genomic annotations and NCI60 processed data at either RNA or protein level
load("data.rda")


# retrieve IPI-HGNC relations for genes profiled at the protein level 
symbol.start = as.numeric(sapply(prot$Fasta.headers, function(x)(regexpr("Gene_Symbol",x)[1])))
ss = substr(prot$Fasta.headers[symbol.start != -1], start=symbol.start[symbol.start != -1], stop=1000000L)
tmp = gsub("Gene_Symbol=","",as.vector(sapply(ss, function(x)(strsplit(x," ")[[1]][1]))))
symbol = as.vector(sapply(tmp, function(x)(strsplit(x,";")[[1]][1])))
zrs = rep(0,length(symbol))
ipiHgncProt = data.frame(IPI=zrs,HGNC.symbol=zrs)
ipiHgncProt[,1] = prot$Accession[symbol.start != -1]
ipiHgncProt[,2] = symbol
# retrieve IPI-HGNC relations for genes profiled at the mRNA level
ipiHgncRna = unique(rna[,2:3])
colnames(ipiHgncRna) = c("IPI","HGNC.symbol")
# retrieve IPI-HGNC relations for genes profiled at either level
ipiHgnc = unique(rbind(ipiHgncProt,ipiHgncRna))


# retrieve genomic annotation for NCI60 profiled genes
l = dim(ipiHgnc)[1]
annotation = data.frame(ID = character(1), IPI = character(l), HGNC.symbol = character(l), Ensembl.Gene.ID = character(l), 
 Chromosome.Name =  character(l), Gene.Start..bp. = integer(l), Gene.End..bp. = integer(l), Strand = character(l))
annotation$ID = paste("ID",1:l,sep="")
annotation$IPI = ipiHgnc[,1] 
annotation$HGNC.symbol = ipiHgnc[,2]
annotation$Ensembl.Gene.ID = biomart$Ensembl.Gene.ID[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Chromosome.Name = biomart$Chromosome.Name[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Gene.Start..bp. = biomart$Gene.Start..bp.[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Gene.End..bp. = biomart$Gene.End..bp.[match(ipiHgnc[,2],biomart$HGNC.symbol)]
annotation$Strand = biomart$Strand[match(ipiHgnc[,2],biomart$HGNC.symbol)]
# remove genes with incomplete genomic annotation
ind = apply(annotation[,4:8], 1, function(x) all(is.na(x)))
annotation = annotation[ !ind, ]

save(annotation,file="annotation.rda")









