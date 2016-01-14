### ==============================================
### MultiAssayExperiment object
### ----------------------------------------------

#' An integrative MultiAssay class for experiment data
#' 
#' @slot Elist A \code{\linkS4class{Elist}} class object for each assay dataset. 
#' @slot pData A \code{DataFrame} of all clinical data available across experiments.
#' @slot sampleMap A \code{DataFrame} of translatable identifiers of samples and participants.
#' @slot metadata Additional data describing the \code{\link{MultiAssayExperiment}} object. 
#' @slot drops A metadata \code{list} of dropped information.
#' @exportClass MultiAssayExperiment
#' @include Elist-class.R
setClass("MultiAssayExperiment",
		 slots = list(
					  Elist="Elist",
					  pData = "DataFrame",
					  sampleMap = "DataFrame",
					  metadata = "ANY",
					  drops = "list"
					  )
		 )

##
## Validity ---------------------------------
##

.uniqueSortIdentical <- function(charvec1, charvec2){
  listInput <- list(charvec1, charvec2)
  listInput <- lapply(listInput, function(x) sort(unique(x)))
  return(identical(listInput[[1]], listInput[[2]]))
}

.allIn <- function(charvec1, charvec2){
  return(all(charvec2 %in% charvec1))
}

## Prepare for unit testing - Travis CI
## sampleMap is a DataFrame with unique sampleNames across assay
.checkSampleMapNames <- function(object){
	errors <- character()
	if(!(.allIn(rownames(pData(object)), slot(sampleMap(object)[, "master"], "values")))){
		msg <- "All samples in the sampleMap must be in the pData"
		errors <- c(errors, msg)
	}
	if(length(errors) == 0L) NULL else errors
}

.uniqueNamesInAssays <- function(object){
  errors <- character()
	SampMap <- sampleMap(object)
	lcheckdups <- S4Vectors::split(SampMap[["assay"]], SampMap[, "assayname"])
	logchecks <- any(vapply(lcheckdups, function(x) any(duplicated(x)), logical(1)))
	if(logchecks){
		msg <- "All sample identifiers in the assays must be unique"
		errors <- c(errors, msg)
	}
	if(length(errors) == 0L) NULL else errors
}

## Experiment list must be the same length as the unique sampleMap assaynames 
.checkElist2 <- function(object){
	errors <- character()
	assaynames <- unique(sampleMap(object)[, "assayname"])
	if(length(Elist(object)) != length(assaynames)){
		msg <- "Elist must be the same length as the sampleMap assaynames"
		errors <- c(errors, msg)
	}
	if(!all(assaynames %in% names(object))){
	  msg <- "Experiment/Assay names in both the Elist and the sampleMap must match"
	  errors <- c(errors, msg)
	}
	if(length(errors) == 0L) NULL else errors
}

## All sample names in the Elist must be in the sampleMap
.checkSampleNames <- function(object){
	if(!identical(sort(unname(unlist(colnames(object)))), 
				  sort(sampleMap(object)[, "assay"]))){
		return("samples in the Elist are not the same as samples in the sampleMap")
	}
	NULL
}

## All names must match between Elist and sampleMap
.checkNames <- function(object){
	if(!all(names(Elist(object)) %in% unique(sampleMap(object)[, "assayname"]))){
		return("Experiment names must match in both Elist and sampleMap")
	}
	NULL
}

.validMultiAssayExperiment <- function(object){
	if(length(Elist(object)) != 0L){
		c(.checkNames(object),	
		  .checkElist2(object), 
		  .checkSampleMapNames(object),
		  .uniqueNamesInAssays(object),
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
		  c_elist <- class(Elist(object))
		  c_mp <- class(pData(object))
		  c_sm <- class(sampleMap(object))
		  cat(sprintf('A "%s"', o_class),
			  "object of", o_len, "listed\n", ifelse(o_len == 1L, "experiment", "experiments"), 
			  "with", ifelse(all(o_names == "none"), "no user-defined names",
							 ifelse(length(o_names) == 1L, "a user-defined name", "user-defined names")),
			  ifelse(length(o_len) == 0L, "or", "and"),
			  ifelse(length(o_len) == 0L, "classes.",
					 ifelse(length(classes) == 1L, "respective class.", "respective classes.")),
			  "\n Containing an ") 
		  show(Elist(object))
		  cat("To access slots use: \n Elist() - to obtain the", sprintf('"%s"', c_elist), 
			  "of experiment instances", 
			  "\n pData() - for the phenotype", sprintf('"%s"', c_mp), 
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
#' @return A \code{\linkS4class{Elist}} object.
#' @exportMethod Elist
#' @describeIn MultiAssayExperiment Access Elist class from MultiAssayExperiment
setMethod("Elist", "MultiAssayExperiment", function(x)
		  getElement(x, "Elist"))

#' Generic Accessor Functions
#' 
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @return A \code{DataFrame} object
#' @exportMethod pData
setGeneric("pData", function(x) standardGeneric("pData"))
#' @describeIn MultiAssayExperiment Access pData slot from MultiAssayExperiment
setMethod("pData", "MultiAssayExperiment", function(x)
		  getElement(x, "pData"))

#' Generic Acessor Functions 
#' 
#' @return Any type of object describing the metadata
#' @describeIn MultiAssayExperiment Access metadata slot from MultiAssayExperiment
#' @exportMethod metadata
setMethod("metadata", "MultiAssayExperiment", function(x)
		  getElement(x, "metadata"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' Length of Experiment List
#' @return An \code{integer} 
#' @exportMethod length
#' @describeIn MultiAssayExperiment Get the length of Elist 
setMethod("length", "MultiAssayExperiment", function(x)
  length(getElement(x, "Elist"))
)

#' @return A character vector of experiment names
#' @exportMethod names
#' @describeIn MultiAssayExperiment Get the names of the Elist
setMethod("names", "MultiAssayExperiment", function(x)
  names(getElement(x, "Elist"))
)

#' @exportMethod
setGeneric("sampleMap<-", function(x, value) standardGeneric("sampleMap<-"))
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "DataFrame"),
                 function(x, value) {
                   slot(x, "sampleMap") <- value
                   return(x)
                 })

#' @exportMethod
setGeneric("Elist<-", function(x, value) standardGeneric("Elist<-"))
setReplaceMethod("Elist", c("MultiAssayExperiment", "Elist"),
                 function(x, value) {
                   slot(x, "Elist") <- value
                   return(x)
                 })
