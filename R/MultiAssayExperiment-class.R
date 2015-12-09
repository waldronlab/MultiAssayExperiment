### ==============================================
### MultiAssayExperiment object
### ----------------------------------------------

#' An integrative MultiAssay class for experiment data
#' 
#' @slot Elist A \code{\linkS4class{Elist}} class object for each assay dataset. 
#' @slot masterPheno A \code{data.frame} of all clinical data available across experiments.
#' @slot sampleMap A \code{data.frame} of translatable identifiers of samples and participants.
#' @slot metadata Additional data describing the \code{\link{MultiAssayExperiment}} object. 
#' @slot drops A metadata \code{list} of dropped information.
#' @exportClass MultiAssayExperiment
#' @include Elist-class.R
setClass("MultiAssayExperiment",
		 slots = list(
					  Elist="Elist",
					  masterPheno = "DataFrame",
					  sampleMap = "DataFrame",
					  metadata = "ANY",
					  drops = "list"
					  )
		 )

##
## Validity ---------------------------------
##

## masterPheno should always be a data.frame
.checkMasterPheno <- function(object){
	errors <- character()
	if(!is(object@masterPheno, "DataFrame")){
		return("masterPheno should be a DataFrame of metadata for all samples!")
	}
	NULL
}

## sampleMap is a DataFrame with unique sampleNames across assay
.checkSampleMap <- function(object){
	errors <- character()
	if(!all(unique(object@sampleMap[, "master"]) %in% rownames(object@masterPheno))){
		msg <- paste("All samples in the sampleMap must be in the masterPheno!")
		errors <- c(errors, msg)
	}
	if(!length(object@Elist) == length(names(S4Vectors::split(object@sampleMap, object@sampleMap[, "assayname"])))){
		msg <- paste("assaynames must be of the same length as the Elist!")
		errors <- c(errors, msg)
	}
	lcheckdups <- S4Vectors::split(object@sampleMap[["assay"]], object@sampleMap$assayname)
	logchecks <- any(vapply(lcheckdups, function(x) any(duplicated(x)), logical(1)))
	if(logchecks){
		msg <- paste("All sample identifiers in the assays must be unique!")
		errors <- c(errors, msg)
	}
	if(length(errors) == 0L) NULL else errors 
}

## Experiment list must be the same length as the unique sampleMap assaynames 
.checkElist2 <- function(object){
	errors <- character()
	assaynames <- unique(object@sampleMap[, "assayname"])
	if(length(object@Elist) != length(assaynames)){
		msg <- paste("Elist must be the same length as the sampleMap assaynames!")
		errors <- c(errors, msg)
	}
	if(!all(assaynames %in% names(object))){
	  msg <- paste("Experiment/Assay names in both the Elist and the sampleMap must match!")
	  errors <- c(errors, msg)
	}
	if(length(errors) == 0L) NULL else errors
}

## All sample names in the Elist must be in the sampleMap
.checkSampleNames <- function(object){
	if(!all.equal(unname(unlist(samples(object))), 
				  object@sampleMap[, "assay"])){
		return("samples in the Elist are not the same as samples in the sampleMap!")
	}
	NULL
}

## All names must match between Elist and sampleMap
.checkNames <- function(object){
	if(!all(names(object@Elist) %in% unique(object@sampleMap[, "assayname"]))){
		return("Experiment names must match in both Elist and sampleMap!")
	}
	NULL
}

.validMultiAssayExperiment <- function(object){
	if(length(object@Elist) != 0L){
		c(.checkMasterPheno(object), 
		  .checkNames(object),	
		  .checkElist2(object), 
		  .checkSampleMap(object),
		  .checkSampleNames(object)
		  )
}}

S4Vectors::setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)


#' Show method for \code{\linkS4class{MultiAssayExperiment}} class
#' 
#' @param object A \code{\link{MultiAssayExperiment}} object.
#' @return Returns a summary of contents for the \code{\link{MultiAssayExperiment}} class. 
#' @exportMethod show
setMethod("show", "MultiAssayExperiment", function(object){
		  o_class <- class(object)
		  o_len <- length(object)
		  o_names <- names(object)
		  if(length(o_names)==0L){ o_names <- "none" }
		  classes <- vapply(Elist(object), class, character(1))
		  c_elist <- class(object@Elist)
		  c_mp <- class(object@masterPheno)
		  c_sm <- class(object@sampleMap)
		  c_md <- class(object@metadata)
		  cat(sprintf('A "%s"', o_class),
			  "object of", o_len, "listed\n", ifelse(o_len == 1L, "experiment", "experiments"), 
			  "with", ifelse(all(o_names == "none"), "no user-defined names",
							 ifelse(length(o_names) == 1L, "a user-defined name", "user-defined names")),
			  ifelse(length(o_len) == 0L, "or", "and"),
			  ifelse(length(o_len) == 0L, "classes.",
					 ifelse(length(classes) == 1L, "respective class.", "respective classes.")),
			  "\n Containing an ") 
		  show(object@Elist)
		  cat("To access slots use: \n Elist() - to obtain the", sprintf('"%s"', c_elist), 
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
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object.
#' @return A \code{list} object. 
#' @exportMethod sampleMap
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))
#' @describeIn MultiAssayExperiment Access sampleMap slot from MultiAssayExperiment
setMethod("sampleMap", "MultiAssayExperiment", function(x)
		  getElement(x, "sampleMap"))

#' Generic Accessor Functions
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object.
#' @return A \code{\link[S4Vectors]{SimpleList-class}} object.
#' @exportMethod Elist
#' @describeIn MultiAssayExperiment Access Elist class from MultiAssayExperiment
setMethod("Elist", "MultiAssayExperiment", function(x)
		  getElement(x, "Elist"))

#' Generic Accessor Functions
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return A \code{data.frame} object
#' @exportMethod masterPheno
setGeneric("masterPheno", function(x) standardGeneric("masterPheno"))
#' @describeIn MultiAssayExperiment Access masterPheno slot from MultiAssayExperiment
setMethod("masterPheno", "MultiAssayExperiment", function(x)
		  getElement(x, "masterPheno"))

#' Generic Acessor Functions 
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return Any type of object describing the metadata
#' @describeIn MultiAssayExperiment Access metadata slot from MultiAssayExperiment
#' @exportMethod metadata
setMethod("metadata", "MultiAssayExperiment", function(x)
		  getElement(x, "metadata"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' Length of Experiment List
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return An \code{integer} 
#' @exportMethod length
#' @describeIn MultiAssayExperiment Get the length of Elist 
setMethod("length", "MultiAssayExperiment", function(x)
  length(getElement(x, "Elist"))
)

#' Names of Experiments 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn MultiAssayExperiment Get the names of the Elist
setMethod("names", "MultiAssayExperiment", function(x)
  names(getElement(x, "Elist"))
)


