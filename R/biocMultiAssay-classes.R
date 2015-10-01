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

#' Feature extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
#' @exportMethod "featExtractor"
#' @export "featExtractor"
setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
setMethod("featExtractor", "ExpressionSet", function(x) affy::featureNames(x))
setMethod("featExtractor", "SummarizedExperiment", function(x) rownames(x))
setMethod("featExtractor", "matrix", function(x) rownames(x))
setMethod("featExtractor", "GRangesList", function(x) ranges(x))

#' Sample extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns an object of the same class  
#' @export "sampleExtractor"
#' @exportMethod "sampleExtractor"
setGeneric("sampleExtractor", function(x) standardGeneric("sampleExtractor"))
setMethod("sampleExtractor", "ExpressionSet", function(x) affy::sampleNames(x)) 
setMethod("sampleExtractor", "SummarizedExperiment", function(x) colnames(x))
setMethod("sampleExtractor", "matrix", function(x) colnames(x))
setMethod("sampleExtractor", "GRangesList", function(x) names(x)) 

#' Subset by Sample method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{"numeric"} or \code{"character"} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @export "subsetSample"
#' @exportMethod "subsetSample"
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
setMethod("subsetSample", "ExpressionSet", function(x, j) x[, j])
setMethod("subsetSample", "SummarizedExperiment", function(x, j) x[, j])
setMethod("subsetSample", "GRangesList", function(x, j) x[j]) 

#' Subset by Feature method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{"numeric"} or \code{"character"} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returnss a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @export "subsetFeature"
#' @exportMethod "subsetFeature"
setGeneric("subsetFeature", function(x, j, ...) standardGeneric("subsetFeature"))
setMethod("subsetFeature", "matrix", function(x, j) x[j, , drop = FALSE])
setMethod("subsetFeature", "ExpressionSet", function(x, j) x[j, ])
setMethod("subsetFeature", "SummarizedExperiment", function(x, j) x[j, ])
# setMethod("subsetFeature", "GRangesList", function(x, j) )

