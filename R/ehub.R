
setClass("expt", representation(tag="character", serType="character",
  assayPath="character", sampleDataPath="character"))

# intent is to refer to serialization type of experimental data,
# to the path by which the assay data is stored, and to
# the path of the sample data

setClass("eHub", representation(hub="list", metadata="ANY", allids="character",
  masterSampleData="data.frame"))
setMethod("phenoData", "eHub", function(object)
  object@masterSampleData)
## setValidity("eHub", function(object){
##   ## Note - requiring existence of local files too restrictive,
##   ## should we allow off-site files e.g. through AnnotationHub
##   ## or other remote file services?
##    all.files <- sapply(object@hub, function(x) x@assayPath)
##     if(!all(file.exists(all.files))){
##       msg <- paste("The following files are not found:",
##                    all.files[!file.exists(all.files)], collapse=", ")
##     }else{
##       return(TRUE)
##     }
##   }
## )

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

createHub <- function(masterpheno, objlist, drop=FALSE, samplemaps=NULL){
  ## samplemaps will be maps that rename samples in object list to names used in masterpheno.
  if(!is(masterpheno, "data.frame"))
     stop("masterpheno should be a data.frame of metadata for all samples")
  if(!is(objlist, "list"))
      stop("objlist should be a named list of data objects")
  ##-----------------------
  ##TODO: sample names mapping if samplemaps provided
  ##-----------------------
  ##Sample names checking:
  has.pheno <- lapply(objlist, function(x) colnames(x) %in% rownames(masterpheno))
  if(!drop){
      errmsg <- paste("Missing the following number of masterpheno entries for each data type: ",
                      paste(names(objlist), ":", sapply(has.pheno, function(x) sum(!x)), collapse=", "),
                      ". Set drop=TRUE to drop these observations, or add samples to masterpheno.")
      stop(errormsg)
  }else{
      message("Dropping the following samples:")
      for (i in 1:length(objlist)){
          if(all(has.pheno[[i]])) next
          message(paste(names(objlist)[i], ":", collapse=""))
          message(paste(colnames(objlist[[i]])[!has.pheno[[i]]], collapse=" "))
          message("\n ")
          objlist[[i]] <- objlist[[i]][, has.pheno[[i]]]
      }

  }
  exptlist <- lapply(1:length(objlist), function(i) new("expt",
     serType="in-memory", assayPath="", tag=names(objlist)[i]))
  hub <- new("eHub", hub=exptlist, masterSampleData=masterpheno)
  res <- new("loadedHub", basehub=hub, elist=objlist)
}
