
setClass("expt", representation(tag="character", serType="character",
  assayPath="character", sampleDataPath="character"))

# intent is to refer to serialization type of experimental data,
# to the path by which the assay data is stored, and to
# the path of the sample data

setClass("eHub", representation(hub="list", metadata="ANY", allids="character",
  masterSampleData="data.frame"))

setMethod("show", "eHub", function(object) cat("eHub with",
  length(object@hub), "experiments.\n"))

setGeneric("loadHub", function(hub)standardGeneric("loadHub"))
setMethod("loadHub", "eHub", function(hub) {
  obj = lapply(hub@hub, function(x) get(load(x@assayPath)))
  names(obj) = sapply(hub@hub, function(x) x@tag)
  new("loadedHub", basehub=hub, elist=obj)
})

setClass("loadedHub", representation(basehub="eHub", elist="list"))
setMethod("show", "loadedHub", function(object) {
 cat("loadedHub instance with experiments having tags:\n")
 cat(paste(names(object@elist), sep=", "))
 print(sapply(object@elist, dim))
})
# lapply(object@elist, show) })
