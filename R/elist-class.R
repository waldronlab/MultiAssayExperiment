.getErrors <- function(object, my_fun){
	obj_cl <- class(object)
	e_class <- class(try(get(my_fun)(object), silent = TRUE))
	if(e_class == "try-error"){
		msg <-  paste0("class ", obj_cl, " should have a '", my_fun, "' method!")
		return(msg)
	}
	NULL
}

### ==============================================
### elist class
### ----------------------------------------------

#' An integrative container for assay data
#' @inheritParams S4Vectors::SimpleList
#' @exportClass elist
setClass("elist", contains = "SimpleList")

#' Generic Builder and Accessor Function
setGeneric("elist", function(x) standardGeneric("elist"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert to \code{\link[S4Vectors]{SimpleList-class}}
#'
#' @param x A \code{list} object
#' @return An \code{\linkS4class{elist}} class object
#' @exportMethod elist
#' @describeIn elist Convert a list to a SimpleList to an elist
setMethod("elist", "list",
		  function(x) new("elist", S4Vectors::SimpleList(x)))
#' @describeIn elist Convert a SimpleList to an elist
setMethod("elist", "SimpleList",
		  function(x) new("elist", x))

##
## Validity ---------------------------------
##

.checkElist <- function(object){
	if(length(object) != 0L){
		errors <- character()
		for(i in seq_along(object)){
			samp_err <- .getErrors(object[[i]], "samples")
			if(!is.null(samp_err)){
				errors <- c(errors, paste0("Element [", i, "] of ", samp_err))
			}
			feat_err <- .getErrors(object[[i]], "features")
			if(!is.null(feat_err)){
				errors <- c(errors, paste0("Element [", i, "] of ", feat_err))
			}
			brack_err <- .getErrors(object[[i]], "[")
			if(!is.null(brack_err)){
				errors <- c(errors, paste0("Element [", i, "] of ", brack_err))
			}
		}
		if(length(errors) == 0L){
			NULL
		} else { errors }
	} else { NULL }
}

.validElist <- function(object){
	.checkElist(object)
}

##  Make sure elist is valid before checking all of the sample names 

S4Vectors::setValidity2("elist", .validElist)

#' Show method for \code{\linkS4class{elist}} class
#'
#' @param object A \code{\linkS4class{elist}} object. 
#' @return Returns a summary of contents for the \code{\linkS4class{elist}} class
#' exportMethod show
setMethod("show", "elist", function(object){
		  o_class <- class(object)
		  elem_cl <- vapply(object, class, character(1))
		  o_len <- length(object)
		  o_names <- names(object)
		  sampdim <- vapply(object, FUN = function(obj) { length(samples(obj)) }, FUN.VALUE = integer(1))
		  featdim <- vapply(object, FUN = function(obj) { length(features(obj)) }, FUN.VALUE = integer(1))
		  cat(sprintf('"%s"', o_class), "class object of length", paste0(o_len, ':'),
			  sprintf('\n [%i] %s: "%s" - %s samples, %s features', seq(o_len), o_names, elem_cl, sampdim, featdim), "\n") 
		  })

