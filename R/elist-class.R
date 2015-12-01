.getMethErr <- function(object, my_fun){
	obj_cl <- class(object)
	e_class <- class(try(get(my_fun)(object), silent = TRUE))
	if(e_class == "try-error"){
		msg <-  paste0("class ", obj_cl, " should have a '", my_fun, "' method!")
		return(msg)
	}
	NULL
}

.getNameErr <- function(object){
  obj_cl <- class(object)
  if(obj_cl == "RangedSummarizeExperiment"){
    if(is.null(names(object))){
      msg <- paste0("names in RangedSummarizedExperiment are NULL!")
      return(msg)
    } else { NULL } 
  } else { NULL }
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

.checkMethods <- function(object){
  errors <- character()
  for(i in seq_along(object)){
    samp_err <- .getMethErr(object[[i]], "samples")
    if(!is.null(samp_err)){
      errors <- c(errors, paste0("Element [", i, "] of ", samp_err))
    }
    feat_err <- .getMethErr(object[[i]], "features")
    if(!is.null(feat_err)){
      errors <- c(errors, paste0("Element [", i, "] of ", feat_err))
    }
    brack_err <- .getMethErr(object[[i]], "[")
    if(!is.null(brack_err)){
      errors <- c(errors, paste0("Element [", i, "] of ", brack_err))
    }
  }
  if(length(errors) == 0L){
    NULL
  } else { errors }
}

.checkElistNames <- function(object){
  errors <- character()
  for(i in seq_along(object)){
    name_err <- .getNameErr(object[[i]])
    if(!is.null(name_err)){
      errors <- c(errors, paste0("[", i, "] Element", name_err))
    }
  }
  if(any(duplicated(names(object)))){
    msg <- paste("Non-unique names provided!")
    errors <- c(errors, msg)
  } 
  if(length(errors) == 0L){
    NULL
  } else { errors }
}

.validElist <- function(object){
  if(length(object) != 0L){
  c(.checkMethods(object),
    .checkElistNames(object))
  }
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

