# .getErrors <- function(object, FUN, ...){
# 	ob_cl <- class(object)
# 	err_cl <- class(try(FUN(object), silent = TRUE))
# 	names(err_cl) <- ob_cl
# 	unsup <- Filter(function(x) x == "try-error", err_cl)
# 	if(err_cl == "try-error"){
# 		return(cat(sprintf("class %s must have a %s method!", names(err_cl), match.call(FUN)[2])))
# 	}
# }

### ==============================================
### elist class
### ----------------------------------------------

#' An integrative container for assay data
#' @exportClass elist
setClass("elist", contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert to \code{\link[S4Vectors]{SimpleList-class}}
#' @param x A \code{list} object
#' @return An \code{\linkS4class{elist}} class object
#' @exportMethod elist
#' @describeIn elist
setMethod("elist", "list",
		  function(x) new("elist", S4Vectors::SimpleList(x)))
#' @describeIn elist
setMethod("elist", "SimpleList",
		  function(x) new("elist", x))


##
## Validity ---------------------------------
##

.checkElist <- function(object){
	if(length(object) != 0L){
		errors <- character()
		objcl <- sapply(object, class)
		featclasses <- sapply(lapply(object, FUN = function(explist) {try(features(explist), silent = TRUE)}), class) 
		featerrors <- featclasses == "try-error"
		sampclasses <- sapply(lapply(object, FUN = function(explist) {try(samples(explist), silent = TRUE)}), class) 
		samperrors <- sampclasses == "try-error" 
		brackcl <- sapply(lapply(object, FUN = function(explist) {try(explist[1], silent = TRUE)}), class)
		brackerr <- brackcl == "try-error"
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
		if(any(brackerr)){
			index <- which(brackcl == "try-error")
			unsupport <- objcl[brackerr]
			msgs <- sapply(seq_along(index), function(x, i) { paste0("Element [", x[i], "] of class '", unsupport[i], "' in the elist must have a bracket '[' method!") }, x = index)
			errors <- c(errors, unname(msgs))
		}
		if(length(errors) == 0L) NULL else errors
	}
	NULL
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
		  cat(sprintf('A "%s"', o_class), "object of length", o_len,
			  sprintf('\n [%i] %s: "%s" - %s samples, %s features', seq(o_len), o_names, elem_cl, sampdim, featdim), "\n") 
		  })

