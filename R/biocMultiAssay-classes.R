.assertExperiment <- function(object) {
    if(!is(object, "SerializedExperiment") && !is(object, "LoadedExperiment"))
        stop("'object' needs to be of classes 'LoadedExperiment' or 'SerializedExperiment'")
}

.assertMultiAssayExperiment <- function(object) {
    if(!is(object, "MultiAssayExperiment"))
        stop("'object' needs to be of class 'MultiAssayExperiment'")
}

.assertScalar <- function(x) {
    if(!is.vector(x) && length(x) == 1 && !is.list(x))
        stop("'x' needs to be a scalar")
}

getExperiments <- function(object) {
    .assertMultiAssayExperiment(object)
    object@elist
}

#' An S4 class for storing experiment data
#' 
#' @slot tag Informal labels for constituents
#' @slot serType Data file type
#' @slot assayPath Assay data location path 
#' @slot sampleDataPath
setClass("expt", representation(tag = "character", serType = "character", 
	assayPath = "character", sampleDataPath = "character")) 

#' An S4 class for storing multiple experiment objects
#' 
#' @slot hub A list of experiment objects 
#' @slot metadata Slot for storing data
#' @slot allids A character vector for patient IDs
#' @slot masterSampleData A data frame for storing clinical data
setClass("eHub", representation(hub = "list", metadata = "ANY", allids = "character", 
	masterSampleData = "data.frame"))

#' Get phenotype data method for eHub class
#' 
#' @param object An \code{\linkS4class{eHub}} class object
setMethod("phenoData", "eHub", function(object) object@masterSampleData)
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

#' Show method for eHub class
#' 
#' @param object An \code{\linkS4class{eHub}} class object
#' @return Returns a list of contents for the eHub class
#' 
setMethod("show", "eHub", function(object){
  cat("eHub with", length(object@hub),
       "experiments.  User-defined tags:\n")
  tags = sapply(object@hub, slot, "tag")
  for (i in 1:length(tags)) {
    cat("\t", tags[i], "\n")
  }
  pd = phenoData(object)
  cat("Sample level data is ", nrow(pd), " x ", ncol(pd), ".\n", sep="")
})

#' Load method for eHub to MultiAssayExperiment class
#' 
#' @param hub An \code{\linkS4class{eHub}} class object
#' @rdname eHub
#' @aliases loadHub, eHub, eHub-method
#' @return Returns a \code{\linkS4class{MultiAssayExperiment}} class object 
setGeneric("loadHub", function(hub)standardGeneric("loadHub"))
setMethod("loadHub", "eHub", function(hub) {
  obj = lapply(hub@hub, function(x) get(load(x@assayPath)))
  names(obj) = sapply(hub@hub, function(x) x@tag)
  new("MultiAssayExperiment", basehub=hub, elist=obj)
})

#' Feature extractor for eSet and SummarizedExperiment
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}} or \code{\linkS4class{SummarizedExperiment}} class object
#' @return Returns either rownames or featureNames
setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
setMethod("featExtractor", "ExpressionSet", function(x) featureNames(x))
setMethod("featExtractor", "SummarizedExperiment", function(x) rownames(x))

#' Show method for MultiAssayExperiment class
#' 
#' @param object A \code{\linkS4class{MultiAssayExperiment}} 
#' @return Returns a list of contents for the MultiAssayExperiment
setClass("MultiAssayExperiment", representation(basehub="eHub", elist="list", sampleData = "DataFrame"))
setMethod("show", "MultiAssayExperiment", function(object) {
 dimmat = t(sapply(object@elist, dim))
 colnames(dimmat) = c("Features", "Samples") # dim for eSet nicer than for SE!
 featExemplars = lapply(object@elist, function(x) head(featExtractor(x),3))
 featExemplars = sapply(featExemplars, paste, collapse=", ")
 featExemplars = substr(featExemplars, 1, 25)
 featExemplars = paste(featExemplars, "...")
 dimmat = data.frame(dimmat, feats.=featExemplars)
 print(dimmat)
})

setMethod("getTag", "MultiAssayExperiment", function(object, i) {
	  if(missing(i)){
	      sapply(object@basehub@hub, FUN = function(x) {getElement(x, "tag")})
	  } else { object@basehub@hub[[i]]@tag }
})

#' Subset method for MultiAssayExperiment class
#' 
#' @return Returns a subset of the \code{\linkS4class{MultiAssayExperiment}} object
setMethod("subset", "MultiAssayExperiment", function(object, samples=NULL, exps = NULL, drop = FALSE) {
    .assertMultiAssayExperiment(object)
if(!is.null(samples)){
    if(is.numeric(samples)) {
        samples <- sampleNames(object@elist)[samples]
    } else { 
    object@elist <- lapply(object@elist, function(oo) {
        jj <- samples[samples %in% sampleNames(oo)]
        oo <- oo[,jj,drop = drop]
        oo
    }) 
}
    object@basehub@masterSampleData <- object@basehub@masterSampleData[samples,]
    object
}
if(!is.null(exps)){
    if(is.character(exps)){
	exps <- match(exps, sapply(object@basehub@hub,FUN = function(x) { getElement(x, "tag") }))
    } 
}
new("MultiAssayExperiment", basehub = getExperiments(object)[oo], elist = getExperiments(object)[exps], sampleData = object@basehub@masterSampleData)
})

