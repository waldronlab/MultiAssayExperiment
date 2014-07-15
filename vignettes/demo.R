
ov = dir(system.file("extdata/tcga_ov", package="biocMultiAssay"), full=TRUE)
tags = c("ov RNA-seq", "ov agilent", "ov mirna", "ov affy", "ov CNV gistic",
  "ov methy 450k")
elist = lapply(1:length(ov), function(x) new("expt", 
     serType="RData", assayPath=ov[x], tag=tags[x], sampleDataPath=ov[x]))
ovhub = new("eHub", hub=elist)
lo = loadHub(ovhub)
ovhub@allids = as.character(unique(lapply(lo, sampleNames)))

