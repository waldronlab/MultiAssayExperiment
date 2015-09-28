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

.checkMap <- function(exptChunk, masterPheno){
	allphenos <- all(exptChunk[,1] %in% rownames(masterPheno))
	uniqss <- all(!duplicated(exptChunk[,2]))
    return(allphenos & uniqss)
}

##
## Validity ---------------------------------
## masterPheno should always be a data.frame
.checkMasterPheno <- function(object){
	if(!is(object@masterPheno, "data.frame")){
		return("masterPheno should be a data frame of metadata for all samples!")
		}
	NULL
}

## SampleMap should be a list of 2 column data.frames
.checkSampleMap <- function(object){
	errors <- character()
	if(!all(sapply(object@sampleMap, is.data.frame))){ 
		msg <- paste("sampleMap must be a list of data.frames!")
		errors <- c(errors, msg)
	}
	if(!all(sapply(object@sampleMap, length)==2)){
		msg <- paste("All data.frames in sampleMap must be of length 2!")
		errors <- c(errors, msg)
	}
	if(!all(sapply(object@sampleMap, .checkMap, object@masterPheno))){
		msg <- paste("sampleMap is not passing all checks!")
		errors <- c(errors, msg)
	}
	if(length(errors) == 0) NULL else errors
}

## Experiment list must be the same length as the sampleMaps list.
.checkElist <- function(object){
	if(length(object@elist) != length(object@sampleMap)){
	return("elist must be the same length as the sampleMap!")
	}
	NULL
}

.getIDs <- function(masterPheno, j){
	return(rownames(masterPheno[j, ]))
}

.getSubsetLogical <- function(object, ids){
	sampCols <- list()
	for(i in seq_along(object@sampleMap)){
		sampCols[i] <- object@sampleMap[[i]][,1] %in% ids
	}
	return(sampCols)
}


#' Feature extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
setMethod("featExtractor", "ExpressionSet", function(x) featureNames(x))
setMethod("featExtractor", "SummarizedExperiment", function(x) rownames(x))
setMethod("featExtractor", "matrix", function(x) rownames(x))
setMethod("featExtractor", "GRangesList", function(x) ranges(x))

#' Sample extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns an object of the same class  
setGeneric("sampleExtractor", function(x) standardGeneric("sampleExtractor"))
setMethod("sampleExtractor", "ExpressionSet", function(x) sampleNames(x)) 
setMethod("sampleExtractor", "SummarizedExperiment", function(x) colnames(x))
setMethod("sampleExtractor", "matrix", function(x) colnames(x))
setMethod("sampleExtractor", "GRangesList", function(x) names(x)) 


#' Subset by Sample method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
setMethod("subsetSample", "ExpressionSet", function(x, j) x[, j])
setMethod("subsetSample", "SummarizedExperiment", function(x, j) x[, j])
setMethod("subsetSample", "GRangesList", function(x, j) x[j]) 

#' Subset by Feature method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returnss a subsetted \code{\linkS4class{MultiAssayExperiment}} object
setGeneric("subsetFeature", function(x, j, ...) standardGeneric("subsetFeature"))
setMethod("subsetFeature", "matrix", function(x, j) x[j, , drop = FALSE])
setMethod("subsetFeature", "ExpressionSet", function(x, j) x[j, ])
setMethod("subsetFeature", "SummarizedExperiment", function(x, j) x[j, ])
# setMethod("subsetFeature", "GRangesList", function(x, j) )


