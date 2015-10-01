setClass("expt", representation(tag="character", serType="character",
  assayPath="character", sampleDataPath="character"))

# intent is to refer to serialization type of experimental data,
# to the path by which the assay data is stored, and to
# the path of the sample data

setClass("eHub", representation(hub="list", metadata="ANY", 
  masterSampleData="data.frame"))
setGeneric("phenoData", function(object) standardGeneric("phenoData"))
setMethod("phenoData", "eHub", function(object)
  object@masterSampleData)
# setValidity("eHub", function(object){
#   ## Note - requiring existence of local files too restrictive,
#   ## should we allow off-site files e.g. through AnnotationHub
#   ## or other remote file services?
#    all.files <- sapply(object@hub, function(x) x@assayPath)
#     if(!all(file.exists(all.files))){
#       msg <- paste("The following files are not found:",
#                    all.files[!file.exists(all.files)], collapse=", ")
#     }else{
#       return(TRUE)
#     }
#   }
# )

setMethod("show", "eHub", function(object) {
  cat("eHub with", length(object@hub),
       "experiments.  User-defined tags:\n")
  tags = sapply(object@hub, slot, "tag")
  for (i in 1:length(tags)) {
    cat("\t", tags[i], "\n")
  }
  pd = phenoData(object)
  cat("Sample level data is ", nrow(pd), " x ", ncol(pd), ".\n", sep="")
})

setGeneric("loadHub", function(hub)standardGeneric("loadHub"))
setMethod("loadHub", "eHub", function(hub) {
  obj = lapply(hub@hub, function(x) get(load(x@assayPath)))
  names(obj) = sapply(hub@hub, function(x) x@tag)
  new("MultiAssayExperiment", basehub=hub, elist=obj)
})

