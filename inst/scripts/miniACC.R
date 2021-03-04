## The following is the code used to create this mini dataset from the full ACC dataset.
## The full ACC MultiAssayExperiment was created by the pipeline at
## https://github.com/waldronlab/MultiAssayExperiment.TCGA

library(curatedTCGAData)
library(MultiAssayExperiment)
library(RaggedExperiment)
library(Biobase)

acc <- curatedTCGAData("ACC", "*", dry.run = FALSE)
proturl <- paste0("http://genomeportal.stanford.edu/",
    "pan-tcga/target_selection_send_data",
    "?filename=Allprotein.txt")
protmap <- read.csv(proturl, as.is = TRUE)

RPPAgenes <- Filter(nchar, protmap$Genes)
RPPAgenes <- unlist(strsplit(RPPAgenes, ","))
RPPAgenes <- unique(RPPAgenes)

miniACC <-
    acc[RPPAgenes, , c("RNASeq2GeneNorm", "gistict", "RPPAArray", "Mutations")]
mut <- assay(miniACC[["Mutations"]], i = "Variant_Classification")
mut <- ifelse(is.na(mut) | mut == "Silent", 0, 1)

miniACC[["Mutations"]] <- mut
colData(miniACC) <- colData(miniACC)[, c(1:17, 810:822)]

rpparowData <-
    protmap[match(rownames(miniACC[["RPPAArray"]]), protmap$Genes),]
rpparowData <- AnnotatedDataFrame(rpparowData)
featureData(miniACC[["RPPAArray"]]) <- rpparowData

md <- list(
    title = "Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma",
    PMID = "27165744",
    RPPAfeatureDataURL = paste0("http://genomeportal.stanford.edu/",
        "pan-tcga/show_target_selection_file",
        "?filename=Allprotein.txt"),
    colDataExtrasURL = "http://www.cell.com/cms/attachment/2062093088/2063584534/mmc3.xlsx"
)
metadata(miniACC) <- md

mirna <- acc[["miRNASeqGene"]]
mirna <- mirna[rowSums(assay(mirna) >= 5) >= 5, ]
experimentData(mirna)@abstract <-
    "Note: Rows not having at least 5 counts in at least 5 samples were removed."
miniACC <- c(miniACC,
             list(miRNASeqGene = mirna),
             sampleMap = sampleMap(acc)[sampleMap(acc)$assay == "miRNASeqGene",])

 miniACC[["RNASeq2GeneNorm"]] <-
     as(miniACC[["RNASeq2GeneNorm"]], "SummarizedExperiment")
 miniACC[["RPPAArray"]] <-
     as(miniACC[["RPPAArray"]], "SummarizedExperiment")
 miniACC[["miRNASeqGene"]] <-
     as(miniACC[["miRNASeqGene"]], "SummarizedExperiment")

save(miniACC, file = "data/miniACC.RData", compress = "bzip2")

