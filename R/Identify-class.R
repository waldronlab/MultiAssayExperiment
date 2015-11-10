### ==============================================
### Identify class 
### ==============================================

#' An identifier class used for subsetting
#' 
#' @slot query Any class indicator needed to subset
#' @slot keeps A \code{list} indicating valid matches in each assay
#' @slot drops A \code{list} of excluded information due to subsetting
#' @slot type A \code{character} vector indicating method used to search
setClass("Identify", 
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

.validIdentify <- function(object){
	c(.checkDrops(object))
}

S4Vectors::setValidity2("Identify", .validIdentify)
