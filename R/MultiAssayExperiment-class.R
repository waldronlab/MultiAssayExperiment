### ==============================================
### multiAssayExperiment object
### ----------------------------------------------

#' An integrative MultiAssay class for experiment data
#' 
#' @slot elist A \code{\link[S4Vectors]{SimpleList-class}} of data across different types of assays. 
#' @slot masterPheno A \code{"data.frame"} of all clinical data available across experiments.
#' @slot sampleMap A \code{"list"} of translatable identifiers of samples and participants.
#' @slot metadata Additional data describing the \code{\linkS4class{MultiAssayExperiment}} object. 
#' @exportClass MultiAssayExperiment
setClass("MultiAssayExperiment",
		 representation(
						elist="SimpleList",
						masterPheno = "data.frame",
						sampleMap = "list", 
						metadata = "ANY"
						) 
		 )

##
## Validity ---------------------------------
##
## Unique samples & phenos all present 
.checkMap <- function(exptChunk, masterPheno){
	allphenos <- all(exptChunk[,1] %in% rownames(masterPheno))
	uniqss <- all(!duplicated(exptChunk[,2]))
    return(allphenos & uniqss)
}

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

.checkNames <- function(object){
	if(!identical(names(object@elist), names(object@sampleMap))){
		return("Experiment names must match in both elist and sampleMap!")
	}
	NULL
}


.validMultiAssayExperiment <- function(object){
  c(.checkMasterPheno(object), 
    .checkSampleMap(object),
    .checkElist(object), 
	.checkNames(object))
}

S4Vectors::setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)

#' Show method for MultiAssayExperiment class
#' 
#' param object A \code{\linkS4class{MultiAssayExperiment}} object.
#' return Returns a summary of contents for the \code{\linkS4class{MultiAssayExperiment}} class. 
#' exportMethod "show"
# setMethod("show", "MultiAssayExperiment", function(object){
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


### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###
#' Generic Accessor Functions
#' @param x A \code{\link{MultiAssayExperiment}} object.
#' @return A \code{"list"} object. 
#' @exportMethod sampleMap
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))
#' @describeIn sampleMap
setMethod("sampleMap", "MultiAssayExperiment", function(x)
getElement(x, "sampleMap"))

#' Generic Accessor Functions
#' @param x A \code{\link{MultiAssayExperiment}} object.
#' @return A \code{\link[S4Vectors]{SimpleList-class}} object.
#' @exportMethod elist
setGeneric("elist", function(x) standardGeneric("elist"))
#' @describeIn elist
setMethod("elist", "MultiAssayExperiment", function(x)
getElement(x, "elist"))

#' Generic Accessor Functions
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return A \code{data.frame} object
#' @exportMethod masterPheno
setGeneric("masterPheno", function(x) standardGeneric("masterPheno"))
#' @describeIn masterPheno
setMethod("masterPheno", "MultiAssayExperiment", function(x)
getElement(x, "masterPheno"))


### - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' Length of Experiment List
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return An \code{integer} 
#' @exportMethod length
#' @describeIn MultiAssayExperiment
setMethod("length", "MultiAssayExperiment", 
	function(x) length(getElement(x, "elist"))
)

#' Names of Experiments 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn MultiAssayExperiment
setMethod("names", "MultiAssayExperiment", 
	function(x) names(getElement(x, "elist"))
)
