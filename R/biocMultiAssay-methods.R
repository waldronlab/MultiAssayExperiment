#' @include elist-class.R stage-class.R MultiAssayExperiment-class.R 
NULL

#' Feature extractor method
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
#' @exportMethod features
setGeneric("features", function(x) standardGeneric("features"))
#' @describeIn features Get the featureNames for ExpressionSet
setMethod("features", "ExpressionSet", function(x) Biobase::featureNames(x))
#' @describeIn features Get a summary of rowRanges for RangedSummarizedExperiment
setMethod("features", "RangedSummarizedExperiment", function(x) names(SummarizedExperiment::rowRanges(x)))
#' @describeIn features Get the rownames of a matrix
setMethod("features", "matrix", function(x) rownames(x))
#' @describeIn features Get the names of a GRanges
setMethod("features", "GRanges", function(x) names(x))
#' @describeIn features Get the summary of ranges for a GRangesList
setMethod("features", "GRangesList", function(x) unlist(sapply(seq_along(x), FUN = function(grl, i)
{paste(rep(names(grl)[i], length(grl[[i]])), features(x[[i]]), sep = "///")}, grl = x)))
#' @describeIn features Get all the features for a MultiAssayExperiment
setMethod("features", "MultiAssayExperiment", function(x) lapply(x@elist, features))

#' Sample extractor generic
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns an object of the same class  
#' @exportMethod samples
#' @describeIn samples Get the sampleNames for ExpressionSet
setMethod("samples", "ExpressionSet", function(object) Biobase::sampleNames(object)) 
#' @describeIn samples Get colnames for RangedSummarizedExperiment
setMethod("samples", "RangedSummarizedExperiment", function(object) BiocGenerics::colnames(object))
#' @describeIn samples Get the colnames of a matrix
setMethod("samples", "matrix", function(object) colnames(object))
#' @describeIn samples Get the names of each list element for a GRangesList
setMethod("samples", "GRangesList", function(object) names(object)) 
#' @describeIn samples Get all the samples for a MultiAssayExperiment
setMethod("samples", "MultiAssayExperiment", function(object) lapply(object@elist, samples))

#' Find hits by class type
#' 
#' @param subject Any valid element from the \code{\linkS4class{elist}} class
#' @param query Either a \code{character} vector or \code{\linkS4class{GRanges}} object used to search by name or ranges
#' @param ... Additional arguments to findOverlaps
#' @return Names of matched queries
#' @exportMethod getHits
setGeneric("getHits", function(subject, query, ...) standardGeneric("getHits"))
#' @describeIn getHits Find overlaps and return names
setMethod("getHits", signature("GRanges", "GRanges"), function(subject, query, ...)
 names(subject[findOverlaps(subject, query, ...)@subjectHits]))
#' @describeIn getHits Iteratively find overlaps and return names
setMethod("getHits", signature("GRangesList", "GRanges"), function(subject, query, ...)
 	lapply(subject, function(grel) { names(grel[findOverlaps(grel, query, ...)]) } ))
#' @describeIn getHits Find overlaps and return names for RangedSummarizedExperiment
setMethod("getHits", signature("RangedSummarizedExperiment", "GRanges"), function(subject, query, ...)
  names(subject[findOverlaps(rowRanges(subject), query, ...)@subjectHits]))
setMethod("getHits", signature("character", "GRanges"), function(subject, query, ...) character(0L) )
#' @describeIn getHits Find matching features in ExpressionSet
setMethod("getHits", signature("ExpressionSet", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
#' @describeIn getHits Find matching features in matrix
setMethod("getHits", signature("matrix", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
#' @describeIn getHits Find all matching features by character
setMethod("getHits", signature("MultiAssayExperiment", "character"), function(subject, query, ...)
  lapply(subject, FUN = function(elem) { getHits(elem, query, ...) }))
#' @describeIn getHits Find all matching features by GRanges
setMethod("getHits", signature("MultiAssayExperiment", "GRanges"), function(subject, query, ...)
  lapply(subject, FUN = function(elem) { getHits(elem, query, ...) }))

#' Subset by Sample generic 
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{numeric}, \code{character}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @exportMethod subsetSample
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
#' @describeIn subsetSample Select columns of a matrix
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
#' @describeIn subsetSample Select samples of an ExpressionSet
setMethod("subsetSample", "ExpressionSet", function(x, j) x[, j])
#' @describeIn subsetSample Select column data of a RangedSummarizedExperiment
setMethod("subsetSample", "RangedSummarizedExperiment", function(x, j) x[,j = j])
#' @describeIn subsetSample Select samples for a GRangesList
setMethod("subsetSample", "GRangesList", function(x, j) x[i=j]) 


#' Subset by Feature method
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{numeric}, \code{character}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @export subsetFeature
setGeneric("subsetFeature", function(x, j, ...) standardGeneric("subsetFeature"))
setMethod("subsetFeature", signature("ANY", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature Subset matrix with corresponding class
setMethod("subsetFeature", signature("matrix", "ANY"), function(x, j){
		  if(any(rownames(x) %in% j)){
			  x <- x[intersect(j, rownames(x)), , drop = FALSE]
			  return(x)
		  } else { 
			  return(x[0, ])
		  }
})
setMethod("subsetFeature", signature("matrix", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature Subset ExpressionSet with corresponding class
setMethod("subsetFeature", signature("ExpressionSet", "ANY"), function(x, j){
		  if(any(featureNames(x) %in% j)){
			  x <- x[intersect(j, featureNames(x)), ]
			  return(x)
		  } else {
			  return(x[0,])
		  }
})
setMethod("subsetFeature", signature("ExpressionSet", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature Use subsetByOverlaps with additional and optional arguments
setMethod("subsetFeature", signature("RangedSummarizedExperiment", "GRanges"), function(x, j, ...){
		  return(subsetByOverlaps(x, j, ...))
})
setMethod("subsetFeature", signature("RangedSummarizedExperiment", "ANY"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature Use subsetByOverlaps for all of the GRangesList
setMethod("subsetFeature", signature("GRangesList", "GRanges"), function(x, j, ...){
		  return(endoapply(x, FUN = function(GR) { subsetByOverlaps(GR, j, ...) })) 
})
setMethod("subsetFeature", signature("GRangesList", "ANY"), function(x, j){ 
		  return(endoapply(x, FUN = function(GR) { GR[0, ] }))
})

#' Convert Stage slot "keeps" to Map
#'
#' @param object An \linkS4class{stage} class object
#' @return Returns a data.frame representation of samples
#' @export getMap
setGeneric("getMap", function(object) standardGeneric("getMap"))
#' describeIn getMap Convert map from list to data.frame
setMethod("getMap", "stage", function(object){
		  if(object@type == "samples"){ return(.ldmap(object@keeps)) }
})

#' Names of Experiments 
#' @param x A \code{\link{stage}} class object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn stage Get the names from the kept elements in the elist
setMethod("names", "stage", function(x)
  names(getElement(x, "keeps"))
)

#' @describeIn stage Get the length of the kept slot 
setMethod("length", "stage", function(x)
  length(getElement(x, "keeps"))
  )

#' Generic Accessor Functions
#' @param x A \code{\linkS4class{stage}} class object
#' @return A \code{character} atomic vector
#' @exportMethod type
setGeneric("type", function(object) standardGeneric("type"))
#' @describeIn stage Get the staging type (either by samples, features, assays)
setMethod("type", "stage", function(object)
  getElement(object, "type")
  )

#' @exportMethod query
setGeneric("query", function(object) standardGeneric("query"))
#' @describeIn stage Get the identifiers used for staging
setMethod("query", "stage", function(object)
  getElement(object, "query"))