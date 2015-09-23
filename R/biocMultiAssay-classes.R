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

checkMap <- function(exptChunk){
	allphenos <- ifelse(all(unique(exptChunk[,1]) %in% rownames(masterPheno)), TRUE, FALSE)
	uniqss <- ifelse(all(!duplicated(exptChunk[,2])), TRUE, FALSE)
	if(uniqss & allphenos){
		return(TRUE)
	} else if(!sampes){
		return(FALSE)
	} else if (!phenos) {
		return(FALSE)
	}
}

checkMAE <- function(object){
	errors <- character()
	if(!is(object@masterPheno, "data.frame")){
		msg <- paste("masterPheno should be a data frame of metadata for all samples")
		errors <- c(errors, msg)
	}
	if(!is(object@elist, "list")){
		msg <- paste("objlist should be a named list of data objects!")
		errors <- c(errors, msg)
	}
	if(!is.null(object@sampleMap)){
		if(any(!sapply(object@sampleMap,checkMap))){
			msg <- paste("The sample maps are not passing the required checks!")
			errors <- c(errors, msg)
		}
		if(length(object@elist) != length(object@sampleMap)){
			msg <- paste("objlist must be the same length as the sampleMap") 
		}
	}
	if(length(errors) == 0) TRUE else errors
}

#' An integrative MultiAssay class for experiment data
#' 
#' @slot elist A list of data across different types of assays 
#' @slot masterPheno A data.frame of all clinical data available across experiments
#' @slot sampleMap A list of translatable identifiers of samples and participants
#' @slot metadata Additional data describing the \code{\linkS4class{MultiAssayExperiment}} class 
setClass("MultiAssayExperiment", representation(elist="list", masterPheno = "data.frame",
	sampleMap = "list", metadata = "ANY"), validity = checkMAE)

#' Show method for MultiAssayExperiment class
#' 
#' @param object A \code{\linkS4class{MultiAssayExperiment}} 
#' @return Returns a list of contents for the MultiAssayExperiment
# setMethod("show", "MultiAssayExperiment", function(object) {
# 		  objdim <- lapply(seq_along(object@elist), FUN = function(j, expt) {	
# 						   dd <- matrix(NA, nrow = length(expt), ncol = 2)
# 						   if(any(is(expt[j], "data.frame"), is(expt[j], "matrix"))){
# 							   dimmat <- matrix(c(dim(expt[j])[1], dim(expt[j])[2]), ncol = 2)
# 							   colnames(dimmat) <- c("Features", "Samples")
# 							   dd <- rbind(dd, dimmat)
# 						   }
# 	} , expt = object@elist)
# print(objdim)
# })


#' Feature extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns either rownames or featureNames
setGeneric("featExtractor", function(x) standardGeneric("featExtractor"))
setMethod("featExtractor", "ExpressionSet", function(x) featureNames(x))
setMethod("featExtractor", signature("SummarizedExperiment", "matrix"), function(x) rownames(x))
setMethod("featExtractor", "GRangesList", function(x) ranges(x))

#' Sample extractor for eSet, SummarizedExperiment, matrix, and GRangesList
#' 
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns an object of the same class  
setGeneric("sampleExtractor", function(x) standardGeneric("sampleExtractor"))
setMethod("sampleExtractor", "ExpressionSet", function(x) sampleNames(x)) 
setMethod("sampleExtractor", signature("SummarizedExperiment", "matrix"), function(x) colnames(x))
setMethod("sampleExtractor", "GRangesList", function(x) names(x)) # if from RTCGAToolbox extract


#' Subset by Sample method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
setMethod("subsetSample", "matrix", function(x, j) x[, j, drop = FALSE])
setMethod("subsetSample", signature("ExpressionSet", "SummarizedExperiment"), function(x, j) x[, j])
setMethod("subsetSample", "GRangesList", function(x, j) x[j]) 

#' Subset by Feature method for eSet, SummarizedExperiment, matrix, and GRangesList
#'
#' @param x Either an \code{\linkS4class{ExpressionSet}}, \code{\linkS4class{GRangesList}}, \code{\linkS4class{SummarizedExperiment}} or \code{matrix} class object
#' @return Returnss a subsetted \code{\linkS4class{MultiAssayExperiment}} object
setGeneric("subsetFeature", function(x, j, ...) standardGeneric("subsetFeature"))
setMethod("subsetFeature", "matrix", function(x, j) x[j, , drop = FALSE])
setMethod("subsetFeature", signature("ExpressionSet", "SummarizedExperiment"), function(x, j) x[j, ])
setMethod("subsetFeature", "GRangesList", function(x, j) )


