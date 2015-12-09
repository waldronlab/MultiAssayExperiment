#' @include Elist-class.R Stage-class.R MultiAssayExperiment-class.R 
NULL

#' Feature extractor method
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
#' @exportMethod features
setGeneric("features", function(x) standardGeneric("features"))
#' @describeIn features Return character if character
setMethod("features", "character", function(x) x)
#' @describeIn features Get the featureNames for ExpressionSet
setMethod("features", "ExpressionSet", function(x) Biobase::featureNames(x))
#' @describeIn features Get a summary of rowRanges for RangedSummarizedExperiment
setMethod("features", "RangedSummarizedExperiment", function(x) names(SummarizedExperiment::rowRanges(x)))
#' @describeIn features Get the rownames of a matrix
setMethod("features", "matrix", function(x) rownames(x))
#' @describeIn features Get the names of a GRanges
setMethod("features", "GRanges", function(x) names(x))
#' @describeIn features Get the summary of ranges for a GRangesList
setMethod("features", "GRangesList", function(x) names(unlist(x, use.names = FALSE)))
#' @describeIn features Get all the features for a MultiAssayExperiment
setMethod("features", "MultiAssayExperiment", function(x) lapply(x@Elist, features))

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
setMethod("samples", "MultiAssayExperiment", function(object) lapply(object@Elist, samples))

#' Find hits by class type
#' 
#' @param subject Any valid element from the \code{\linkS4class{Elist}} class
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
  names(unlist(subject, use.names = FALSE))[findOverlaps(subject, query, ...)@subjectHits])
#' @describeIn getHits Find matching features from GRangesList
setMethod("getHits", signature("GRangesList", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
#' @describeIn getHits Find overlaps and return names for RangedSummarizedExperiment
setMethod("getHits", signature("RangedSummarizedExperiment", "GRanges"), function(subject, query, ...)
  names(subject[findOverlaps(rowRanges(subject), query, ...)@subjectHits]))
#' @describeIn getHits Find matching features from RangedSummarizedExperiment
setMethod("getHits", signature("RangedSummarizedExperiment", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
setMethod("getHits", signature("character", "GRanges"), function(subject, query, ...) character(0L) )
#' @describeIn getHits Find matching features in ExpressionSet
setMethod("getHits", signature("ExpressionSet", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
setMethod("getHits", signature("ExpressionSet", "GRanges"), function(subject, query, ...)
  intersect(features(query), features(subject)))
#' @describeIn getHits Find matching features in matrix
setMethod("getHits", signature("matrix", "character"), function(subject, query, ...)
  intersect(query, features(subject)))
setMethod("getHits", signature("matrix", "GRanges"), function(subject, query, ...)
  intersect(features(query), features(subject)))
#' @describeIn getHits Find all matching features by character
setMethod("getHits", signature("MultiAssayExperiment", "character"), function(subject, query, ...)
  lapply(subject@Elist, FUN = function(elem) { getHits(elem, query, ...) }))
#' @describeIn getHits Find all matching features by GRanges
setMethod("getHits", signature("MultiAssayExperiment", "GRanges"), function(subject, query, ...)
  lapply(subject@Elist, FUN = function(elem) { getHits(elem, query, ...) }))

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
#' @param object An \linkS4class{Stage} class object
#' @return Returns a data.frame representation of samples
#' @export getMap
setGeneric("getMap", function(object) standardGeneric("getMap"))
#' describeIn Stage Convert map from list to DataFrame
setMethod("getMap", "Stage", function(object){
		  .convertList(object, "keeps")
})

#' Names of Experiments 
#' @param x A \code{\link{Stage}} class object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn Stage Get the names from the kept elements in the Elist
setMethod("names", "Stage", function(x)
  names(getElement(x, "keeps"))
)

#' @describeIn Stage Get the number of assays from kept slot 
setMethod("length", "Stage", function(x)
  length(getElement(x, "keeps"))
)

#' Generic Accessor Functions
#' @param x A \code{\linkS4class{Stage}} class object
#' @return A \code{character} atomic vector
#' @exportMethod type
setGeneric("type", function(object) standardGeneric("type"))
#' @describeIn Stage Get the staging type (either by samples, features, assays)
setMethod("type", "Stage", function(object)
  getElement(object, "type")
  )

#' @exportMethod query
setGeneric("query", function(object) standardGeneric("query"))
#' @describeIn Stage Get the identifiers used for staging
setMethod("query", "Stage", function(object)
  getElement(object, "query")
)
