setags = c("SE ov expr", "SE ov CN")
sepaths = c("./TCGA_se.rda", "./gistic.cn.rda")
elistSE = lapply(1:length(setags), function(x) new("expt",
   serType="RData", assayPath=sepaths[x], tag=setags[x], sampleDataPath="")) #
ovhub2SE = new("eHub", hub=elistSE, masterSampleData = as(colData(TCGA_se), "data.frame"))
#lov = loadHub(ovhub2)
