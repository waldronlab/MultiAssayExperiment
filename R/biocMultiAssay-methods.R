#' @include Elist-class.R MultiAssayView-class.R RaggedRangedAssay-class.R MultiAssayExperiment-class.R
#' @import BiocGenerics
NULL

#' Rownames extractor method 
#' 
#' \code{rownames} are used to obtain feature names from experiment data
#' 
#' @param x A compatible object with feature data 
#' @return Returns either rownames or featureNames
#' @exportMethod rownames
setMethod("rownames", "character", function(x) x)
setMethod("rownames", "ExpressionSet", function(x) Biobase::featureNames(x))
setMethod("rownames", "RangedSummarizedExperiment", function(x) names(SummarizedExperiment::rowRanges(x)))
setMethod("rownames", "GRanges", function(x) names(x))
setMethod("rownames", "RaggedRangedAssay", function(x) names(unlist(x, use.names = FALSE)))
#' @describeIn MultiAssayExperiment Get all the rownames for a MultiAssayExperiment
setMethod("rownames", "MultiAssayExperiment", function(x) lapply(x@Elist, rownames))

#' Colnames extractor 
#' 
#' \code{colnames} are used to obtain sample names from experiment data
#' 
#' @param x A compatible object with sample data
#' @return Returns a vector of colnames or samplenames
#' @exportMethod colnames
setMethod("colnames", "ExpressionSet", function(x) Biobase::sampleNames(x)) 
setMethod("colnames", "RaggedRangedAssay", function(x) base::names(x)) 
#' @describeIn MultiAssayExperiment Get all the colnames for a MultiAssayExperiment
setMethod("colnames", "MultiAssayExperiment", function(x) lapply(x@Elist, colnames))

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
setMethod("getHits", signature("RaggedRangedAssay", "GRanges"), function(subject, query, ...)
  names(unlist(subject, use.names = FALSE))[findOverlaps(subject, query, ...)@subjectHits])
#' @describeIn getHits Find matching rownames from RaggedRangedAssay
setMethod("getHits", signature("RaggedRangedAssay", "character"), function(subject, query, ...)
  intersect(query, rownames(subject)))
#' @describeIn getHits Find overlaps and return names for RangedSummarizedExperiment
setMethod("getHits", signature("RangedSummarizedExperiment", "GRanges"), function(subject, query, ...)
  names(subject[findOverlaps(rowRanges(subject), query, ...)@subjectHits]))
#' @describeIn getHits Find matching rownames from RangedSummarizedExperiment
setMethod("getHits", signature("RangedSummarizedExperiment", "character"), function(subject, query, ...)
  intersect(query, rownames(subject)))
setMethod("getHits", signature("character", "GRanges"), function(subject, query, ...) character(0L) )
#' @describeIn getHits Find matching rownames in ExpressionSet
setMethod("getHits", signature("ExpressionSet", "character"), function(subject, query, ...)
  intersect(query, rownames(subject)))
setMethod("getHits", signature("ExpressionSet", "GRanges"), function(subject, query, ...)
  intersect(rownames(query), rownames(subject)))
#' @describeIn getHits Find matching rownames in matrix
setMethod("getHits", signature("matrix", "character"), function(subject, query, ...)
  intersect(query, rownames(subject)))
setMethod("getHits", signature("matrix", "GRanges"), function(subject, query, ...)
  intersect(rownames(query), rownames(subject)))
#' @describeIn getHits Find all matching rownames by character
setMethod("getHits", signature("MultiAssayExperiment", "character"), function(subject, query, ...)
  lapply(subject@Elist, FUN = function(elem) { getHits(elem, query, ...) }))
#' @describeIn getHits Find all matching rownames by GRanges
setMethod("getHits", signature("MultiAssayExperiment", "GRanges"), function(subject, query, ...)
  lapply(subject@Elist, FUN = function(elem) { getHits(elem, query, ...) }))

# SEE : 
  # showMethods("findOverlaps", classes = "RaggedRangedAssay", inherited = TRUE) 
#----
## setMethod("getHits", signature("ANY", "GRanges"), function(subject, query, ...){
##   ranged <- all(c("findOverlaps", "subsetByOverlaps") %in% showMethods(class = class(subject)))
##   if(ranged){
##     query <- names(unlist(subject, use.names = FALSE))[findOverlaps(subject, query, ...)@subjectHits])
##     getHits(subject, query)
##   } else {
##     subject[0,]
##   }
## })
## setMethod("getHits", signature("ANY", "character"), function(subject, query, ...){
##   query[query %in% rownames(subject)]
## })

#' Subset by Sample generic 
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{RaggedRangedAssay}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
#' @param j Either a \code{numeric}, \code{character}, or \code{logical} vector class for subsetting
#' @param ... Additional arguments to pass
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' @exportMethod subsetSample
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
#' @describeIn subsetSample Select columns of a matrix
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
#' @describeIn subsetSample Select colnames of an ExpressionSet
setMethod("subsetSample", "ExpressionSet", function(x, j) x[, j])
#' @describeIn subsetSample Select column data of a RangedSummarizedExperiment
setMethod("subsetSample", "RangedSummarizedExperiment", function(x, j) x[,j = j])
#' @describeIn subsetSample Select colnames for a RaggedRangedAssay
setMethod("subsetSample", "RaggedRangedAssay", function(x, j) x[i=j]) 


#' Subset by Feature method
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{RaggedRangedAssay}}, \code{\linkS4class{RangedSummarizedExperiment}} or \code{matrix} class object
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
#' @describeIn subsetFeature Use subsetByOverlaps for all of the RaggedRangedAssay
setMethod("subsetFeature", signature("RaggedRangedAssay", "GRanges"), function(x, j, ...){
		  return(endoapply(x, FUN = function(GR) { subsetByOverlaps(GR, j, ...) })) 
})
setMethod("subsetFeature", signature("RaggedRangedAssay", "ANY"), function(x, j){ 
		  return(endoapply(x, FUN = function(GR) { GR[0, ] }))
})

#' Convert MultiAssayView slot "keeps" to Map
#'
#' @param object An \linkS4class{MultiAssayView} class object
#' @return Returns a data.frame representation of samples
#' @export getMap
setGeneric("getMap", function(object) standardGeneric("getMap"))
#' describeIn MultiAssayView Convert map from list to DataFrame
setMethod("getMap", "MultiAssayView", function(object){
		  .convertList(object, "keeps")
})

#' Names of Experiments 
#' @param x A \code{\link{MultiAssayView}} class object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn MultiAssayView Get the names from the kept elements in the Elist
setMethod("names", "MultiAssayView", function(x)
  names(getElement(x, "keeps"))
)

#' @describeIn MultiAssayView Get the number of assays from kept slot 
setMethod("length", "MultiAssayView", function(x)
  length(getElement(x, "keeps"))
)

#' Generic Accessor Functions
#' @param x A \code{\linkS4class{MultiAssayView}} class object
#' @return A \code{character} atomic vector
#' @exportMethod type
setGeneric("type", function(object) standardGeneric("type"))
#' @describeIn MultiAssayView Get the staging type (either by samples, features, assays)
setMethod("type", "MultiAssayView", function(object)
  getElement(object, "type")
  )

#' @exportMethod query
setGeneric("query", function(object) standardGeneric("query"))
#' @describeIn MultiAssayView Get the identifiers used for staging
setMethod("query", "MultiAssayView", function(object)
  getElement(object, "query")
)
