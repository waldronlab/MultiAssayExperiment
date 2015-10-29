### ==============================================
### multiAssayExperiment object
### ----------------------------------------------

#' An integrative MultiAssay class for experiment data
#' 
#' @slot elist A \code{\link[S4Vectors]{SimpleList-class}} of data across different types of assays. 
#' @slot masterPheno A \code{"data.frame"} of all clinical data available across experiments.
#' @slot sampleMap A \code{"list"} of translatable identifiers of samples and participants.
#' @slot metadata Additional data describing the \code{\link{MultiAssayExperiment}} object. 
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
.checkMap <- function(mappeddf, masterPheno){
	allmapped <- all(mappeddf[,1] %in% rownames(masterPheno))
	allphenos <- all(rownames(masterPheno) %in% mappeddf[,1])
	uniqss <- all(!duplicated(na.omit(mappeddf[,2])))
	return(allmapped & allphenos & uniqss)
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
	if(!all(sapply(object@sampleMap, length) == 2)){
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
	errors <- character()
	if(length(object@elist) != length(object@sampleMap)){
		msg <- paste("elist must be the same length as the sampleMap!")
		errors <- c(errors, msg)
	}
	objcl <- sapply(object@elist, class)
	featclasses <- sapply(lapply(object@elist, FUN = function(explist) {try(features(explist), silent = TRUE)}), class) 
	featerrors <- featclasses == "try-error"
	sampclasses <- sapply(lapply(object@elist, FUN = function(explist) {try(samples(explist), silent = TRUE)}), class) 
	samperrors <- sampclasses == "try-error" 
	if(any(featerrors)){
		index <- which(featclasses == "try-error")
		unsupport <- objcl[featerrors]
		msgs <- sapply(seq_along(index), function(x, i) { paste0("Element [", x[i], "] of class '", unsupport[i], "' in the elist must have a features method!") }, x = index)
		errors <- c(errors, unname(msgs))
	}
	if(any(samperrors)){
		index <- which(sampclasses == "try-error")
		unsupport <- objcl[samperrors]
		msgs <- sapply(seq_along(index), function(x, i) { paste0("Element [", x[i], "] of class '", unsupport[i], "' in the elist must have a samples method!") }, x = index)
		errors <- c(errors, unname(msgs))
	}
	if(length(errors) == 0) NULL else errors
}

.checkSampleNames <- function(object){
	Map(all.equal,
		lapply(object@elist, samples),
		lapply(object@sampleMap, FUN = function(map) {na.omit(map)[,2]}))
}

.checkNames <- function(object){
	if(!all(names(object@elist) %in% names(object@sampleMap))){
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
#' @param object A \code{\link{MultiAssayExperiment}} object.
#' @return Returns a summary of contents for the \code{\link{MultiAssayExperiment}} class. 
#' @exportMethod show
setMethod("show", "MultiAssayExperiment", function(object){
		  o_class <- class(object)
		  o_len <- length(object)
		  o_names <- names(object)
		  classes <- vapply(elist(object), class, character(1))
		  c_elist <- class(object@elist)
		  c_mp <- class(object@masterPheno)
		  c_sm <- class(object@sampleMap)
		  c_md <- class(object@metadata)
		  cat(sprintf('A "%s"', o_class),
			  "object containing", o_len, 
			  "\n listed", ifelse(o_len == 1L, "experiment", "experiments"), 
			  "with", ifelse(length(o_names) == 0L, "no user-defined names",
							 ifelse(length(o_names) == 1L, "a user-defined name", "user-defined names")),
			  ifelse(length(o_names) == 0L, "or", "and"),
			  ifelse(length(o_names) == 0L, "classes.",
					 ifelse(length(o_names) == 1L, "its respective class:", "their respective classes:")), 
			  sprintf('\n [%i] %s - "%s"', seq(length(o_names)), o_names, classes), "\n")
		  cat("To access slots use: \n elist() - to obtain the", sprintf('"%s"', c_elist), 
			  "of experiment instances", 
			  "\n masterPheno() - for the phenotype", sprintf('"%s"', c_mp), 
			  "\n sampleMap() - for the sample availability", sprintf('"%s"', c_sm), 
			  "\n metadata() - for the metadata object of 'ANY' class",
			  "\nSee also: subsetByAssay(), subsetByFeature(), subsetBySample()\n")
	})


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

#' Generic Acessor Functions 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return Any type of object describing the metadata
#' @exportMethod metadata
setGeneric("metadata", function(x) standardGeneric("metadata"))
#' @describeIn metadata
setMethod("metadata", "MultiAssayExperiment", function(x)
		  getElement(x, "metadata"))

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
