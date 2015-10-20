### ==============================================
### Identify class 
### ==============================================

#' An identifier class used for subsetting
#' 
#' @slot logreturn A \code{logical} vector for subsetting
#' @slot drops A \code{list} of excluded information due to subsetting
setClass("Identify", 
		 representation(
						logreturn = "logical",
						drops = "list"
						)
		 )

.checkLogreturn <- function(object){
	if(!is.logical(object@logreturn)){
		return("logreturn must be a logical vector!")
	}
	NULL
}

.checkDrops <- function(object){
	if(length(object@drops) != length(object@logreturn)){
		return("List of dropped information must be the same as the logical vector!")
	}
	NULL
}

.validIdentify <- function(object){
	c(.checkLogreturn(object), 
	  .checkDrops(object))
}

S4Vectors::setValidity2("Identify", .validIdentify)
