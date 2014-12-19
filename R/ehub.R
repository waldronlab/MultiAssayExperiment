
setClass("expt", representation(tag="character", serType="character",
  assayPath="character", sampleDataPath="character"))

# intent is to refer to serialization type of experimental data,
# to the path by which the assay data is stored, and to
# the path of the sample data

setClass("eHub", representation(hub="list", metadata="ANY", allids="character",
  masterSampleData="data.frame"))
setMethod("phenoData", "eHub", function(object)
  object@masterSampleData)

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
  new("loadedHub", basehub=hub, elist=obj)
})

setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
setMethod("featExtractor", "ExpressionSet", function(x) featureNames(x))
setMethod("featExtractor", "SummarizedExperiment", function(x) rownames(x))

setClass("loadedHub", representation(basehub="eHub", elist="list"))
setMethod("show", "loadedHub", function(object) {
 cat("loadedHub instance.\n")
 dimmat = t(sapply(object@elist, dim))
 colnames(dimmat) = c("Features", "Samples") # dim for eSet nicer than for SE!
 featExemplars = lapply(object@elist, function(x) head(featExtractor(x),3))
 featExemplars = sapply(featExemplars, paste, collapse=", ")
 featExemplars = substr(featExemplars, 1, 25)
 featExemplars = paste(featExemplars, "...")
 dimmat = data.frame(dimmat, feats.=featExemplars)
 print(dimmat)
})
