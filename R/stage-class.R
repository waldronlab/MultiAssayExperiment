### ==============================================
### Stage class 
### ==============================================

#' An identifier class used for staging a subset operation
#' 
#' @slot query Any class indicator needed to subset
#' @slot keeps A \code{list} indicating valid matches in each assay
#' @slot drops A \code{list} of excluded information due to subsetting
#' @slot type A \code{character} vector indicating method used to search
#' @exportClass stage
setClass("stage", 
		 representation(query = "ANY",
						keeps = "list",
						drops  = "list", 
						type = "character")
		 )

.checkDrops <- function(object){
	errors <- character()
	if(length(object@keeps) != 0L){
		if(length(object@drops) != length(object@keeps)){
			msg <- paste("List of dropped information must be the same length as the kept information!")
			errors <- c(errors, msg)
		}
	}
	if(length(errors) == 0L) NULL else errors	
}

.validStage <- function(object){
	c(.checkDrops(object))
}

S4Vectors::setValidity2("stage", .validStage)
