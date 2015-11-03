### ==============================================
### Identify class 
### ==============================================

#' An identifier class used for subsetting
#' 
#' @slot inassay A \code{logical} indicator for a match by assay
#' @slot indim A \code{list} indicating valid matches in each assay
#' @slot identifier \code{ANY} object that in used to find matches
#' @slot drops A \code{list} of excluded information due to subsetting
setClass("Identify", 
		 representation(
						 keeps = "list",
						drops  = "list", 
						type = "character"
						)
		 )

.checkDrops <- function(object){
	errors <- character()
	if(length(object@inassay) == 0L){
		if(length(object@drops) != length(object@indim)){
			msg <- paste("List of dropped information must be the same as the logical vector!")
			errors <- c(errors, msg)
		}
	}
	if(length(object@indim) == 0L){
		if(length(object@drops) != length(object@inassay)){
			msg <- paste("List of dropped information must be the same as the logical vector!")
			errors <- c(errors, msg)
		}
	}
	if(length(errors) == 0L) NULL else errors	
}

.validIdentify <- function(object){
	c(.checkDrops(object))
}

S4Vectors::setValidity2("Identify", .validIdentify)
