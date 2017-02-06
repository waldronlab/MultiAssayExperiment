library(MultiAssayExperiment)

if(!file.exists("gbmMAEO.rds")){
  download.file("http://s3.amazonaws.com/multiassayexperiments/gbmMAEO.rds",
                destfile="gbmMAEO.rds")  
}

gbm = readRDS("gbmMAEO.rds")
lit = gbm[, , c("RNASeq2GeneNorm", "gistica", "Methylation", 
                               "miRNAArray")]

png("upsetSamples.png", width=4, height=4, units="in", res=600)
MultiAssayExperiment::upsetSamples(lit, text.scale=1.3, set_size.angles=45, 
                                   mb.ratio=c(.6,.4))
dev.off()
