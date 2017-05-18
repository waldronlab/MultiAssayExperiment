##
## The following is the code used to create this mini dataset from the full ACC dataset.
## The full ACC MultiAssayExperiment was created by the pipeline at
## https://github.com/waldronlab/multiassayexperiment-tcga.
\dontrun{
    ## See www.tinyurl.com/MAEOurls for more pre-built TCGA MultiAssayExperiment objects
    download.file("http://s3.amazonaws.com/multiassayexperiments/accMAEO.rds",
                  destfile = "~/Downloads/accMAEO.rds")
    library(MultiAssayExperiment)
    library(RaggedExperiment) #needed for RaggedExperiment objects by updateObject()
    acc <- readRDS("~/Downloads/accMAEO.rds")
    acc <- updateObject(acc)
    protmap <-
        read.csv(
            "http://genomeportal.stanford.edu/pan-tcga/target_selection_send_data?filename=Allprotein.txt",
            as.is = TRUE
        )

    RPPAgenes <- protmap$Genes
    RPPAgenes <- c(RPPAgenes, unlist(strsplit(RPPAgenes, ",")))
    RPPAgenes <- unique(RPPAgenes)
    RPPAgenes <- RPPAgenes[!RPPAgenes == ""]

    miniACC <-
        acc[RPPAgenes, , c("RNASeq2GeneNorm", "gistict", "RPPAArray", "Mutations")]
    mut <- assays(miniACC[["Mutations"]])[["Variant_Classification"]]
    mut[is.na(mut)] <- 0
    mut[mut == "Silent"] <- 0
    mut[mut != "0"] <- 1
    class(mut) <- "numeric"

    miniACC[["Mutations"]] <- mut
    colData(miniACC) <- colData(miniACC)[, c(1:17, 810:822)]

    library(Biobase)
    rpparowData <-
        protmap[match(rownames(miniACC[["RPPAArray"]]), protmap$Genes),]
    rpparowData <- AnnotatedDataFrame(rpparowData)
    featureData(miniACC[["RPPAArray"]]) <- rpparowData

    md <-
        list(
            title = "Comprehensive Pan-Genomic Characterization of Adrenocortical Carcinoma",
            PMID = "27165744",
            sourceURL = "http://s3.amazonaws.com/multiassayexperiments/accMAEO.rds",
            RPPAfeatureDataURL = "http://genomeportal.stanford.edu/pan-tcga/show_target_selection_file?filename=Allprotein.txt",
            colDataExtrasURL = "http://www.cell.com/cms/attachment/2062093088/2063584534/mmc3.xlsx"
        )
    metadata(miniACC) <- md

    mirna <- acc[["miRNASeqGene"]]
    mirna <-
        mirna[apply(assay(mirna), 1, function(x)
            sum(x >= 5)) >= 5,]
    experimentData(mirna)@abstract <-
        "Note: Rows not having at least 5 counts in at least 5 samples were removed."
    miniACC <- c(miniACC,
                 list(miRNASeqGene = mirna),
                 sampleMap = sampleMap(acc)[sampleMap(acc)$assay == "miRNASeqGene",])

    save(miniACC, file = "miniACC.RData", compress = "bzip2")
}
