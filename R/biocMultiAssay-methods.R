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

#' Feature extractor method
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
#' @exportMethod featExtractor
setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
#' @describeIn featExtractor
setMethod("featExtractor", "ExpressionSet", function(x) affy::featureNames(x))
#' @describeIn featExtractor
setMethod("featExtractor", "RangedSummarizedExperiment", function(x) SummarizedExperiment::rowRanges(x))
#' @describeIn featExtractor
setMethod("featExtractor", "matrix", function(x) rownames(x))
#' @describeIn featExtractor
setMethod("featExtractor", "GRangesList", function(x) x)

#' Sample extractor generic
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @return Returns an object of the same class  
#' @exportMethod sampleExtractor
setGeneric("sampleExtractor", function(x) standardGeneric("sampleExtractor"))
#' @describeIn sampleExtractor 
setMethod("sampleExtractor", "ExpressionSet", function(x) Biobase::sampleNames(x)) 
#' @describeIn sampleExtractor 
setMethod("sampleExtractor", "RangedSummarizedExperiment", function(x) names(x))
#' @describeIn sampleExtractor 
setMethod("sampleExtractor", "matrix", function(x) colnames(x))
#' @describeIn sampleExtractor 
setMethod("sampleExtractor", "GRangesList", function(x) names(x)) 

#' Subset by Sample generic 
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{numeric}, \code{character}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @exportMethod subsetSample
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
#' @describeIn subsetSample
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
#' @describeIn subsetSample
setMethod("subsetSample", "ExpressionSet", function(x, j) x[, j])
#' @describeIn subsetSample
setMethod("subsetSample", "RangedSummarizedExperiment", function(x, j) x[j,])
#' @describeIn subsetSample
setMethod("subsetSample", "GRangesList", function(x, j) x[i=j]) 

#' Subset by Feature method
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{"numeric"}, \code{"character"}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returnss a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @export subsetFeature
setGeneric("subsetFeature", function(x, j, ...) standardGeneric("subsetFeature"))
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("ANY", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("matrix", "ANY"), function(x, j){
		  if(any(rownames(x) %in% j)){
			  x <- x[intersect(j, rownames(x)), , drop = FALSE]
			  return(x)
		  } else { 
			  return(colnames(x))
		  }
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("matrix", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("ExpressionSet", "ANY"), function(x, j){
		  found <- featureNames(x) %in% j
		  if(any(found)){
			  j <- featureNames(x)[found]
			  x <- x[j, ]
			  return(x)
		  } else {
			  return(x[0,])
		  }
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("ExpressionSet", "GRanges"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("RangedSummarizedExperiment", "GRanges"), function(x, j, ...){
		  return(subsetByOverlaps(x, j, ...))
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("RangedSummarizedExperiment", "ANY"), function(x, j){
		  return(x[0, ])
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("GRangesList", "GRanges"), function(x, j, ...){
		  return(lapply(x, FUN = function(GR) { subsetByOverlaps(GR, j, ...) })) 
})
#' @describeIn subsetFeature
setMethod("subsetFeature", signature("GRangesList", "ANY"), function(x, j){ 
		  return(lapply(x, FUN = function(GR) { x[0, ] }))
})

