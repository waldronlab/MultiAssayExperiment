#' @include elist-class.R MultiAssayExperiment-class.R
NULL

.assertMultiAssayExperiment <- function(object) {
    if(!is(object, "MultiAssayExperiment"))
        stop("'object' needs to be of class 'MultiAssayExperiment'")
}

.assertScalar <- function(x) {
    if(!is.vector(x) && length(x) == 1 && !is.list(x))
        stop("'x' needs to be a scalar")
}

#' Feature extractor method
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
#' @exportMethod features
setGeneric("features", function(x) standardGeneric("features"))
#' @describeIn features Get the featureNames for ExpressionSet
setMethod("features", "ExpressionSet", function(x) Biobase::featureNames(x))
#' @describeIn features Get a summary of rowRanges for RangedSummarizedExperiment
setMethod("features", "RangedSummarizedExperiment", function(x) BiocGenerics::unlist(GenomicRanges::rowRanges(x)))
#' @describeIn features Get the rownames of a matrix
setMethod("features", "matrix", function(x) rownames(x))
#' @describeIn features Get the summary of ranges for a GRangesList
setMethod("features", "GRangesList", function(x) BiocGenerics::unlist(x))
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
#' @param j Either a \code{"numeric"}, \code{"character"}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returnss a subsetted \code{\linkS4class{MultiAssayExperiment}} object
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
