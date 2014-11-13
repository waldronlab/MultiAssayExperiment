
# date: 26/03/2014
# aim: save the initial data from NCI60 and basic Ensembl annotation in an R session
# output: data.rda

# genomic annotations from Ensembl Biomart
biomart = read.table("Ensembl.BioMart.v75.human.txt",header=TRUE,as.is=TRUE,sep="\t")
# NCI60 rna processed data
rna = read.table("Table_S8_NCI60_transcriptome.txt",header=TRUE,as.is=TRUE,sep="\t")
rna[,3] = gsub(" \\(includes.*\\)","",rna[,3])
# NCI60 protein processed data
prot = read.table("Table_S3_NCI60_proteome.txt",header=TRUE,as.is=TRUE,sep="\t")
save(biomart,rna,prot,file="data.rda")

