.getMethErr <- function(object, my_fun){
	obj_cl <- class(object)
	e_class <- class(try(get(my_fun)(object), silent = TRUE))
	if(e_class == "try-error"){
		msg <-  paste0("class ", obj_cl, " should have a '", my_fun, "' method!")
		return(msg)
	}
	NULL
}

createNames <- function(object){
  for(i in seq_along(object)){
   names(object[[i]]) <- 1:length(object[[i]]) 
  }
  return(object)
}

.getNameErr <- function(object){
  obj_cl <- class(object)
  if(obj_cl == "RaggedRangedAssay"){
    if(is.null(names(object))){
      msg <- paste("names in", obj_cl, "are NULL!")
      return(msg)
    } else { NULL } 
  } else { NULL }
}

### ==============================================
### Elist class
### ----------------------------------------------

#' An integrative container for assay data
#' @inheritParams S4Vectors::SimpleList
#' @exportClass Elist
setClass("Elist", contains = "SimpleList")

#' Generic Builder and Accessor Function
setGeneric("Elist", function(x) standardGeneric("Elist"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert to \code{\link[S4Vectors]{SimpleList-class}}
#'
#' @param x A \code{list} object
#' @return An \code{\linkS4class{Elist}} class object
#' @exportMethod Elist
#' @describeIn Elist Convert a list to a SimpleList to an Elist
setMethod("Elist", "list",
		  function(x) new("Elist", S4Vectors::SimpleList(x)))
#' @describeIn Elist Convert a SimpleList to an Elist
setMethod("Elist", "SimpleList",
		  function(x) new("Elist", x))
setMethod("Elist", "missing", 
          function(x) new("Elist"))

##
## Validity ---------------------------------
##

.checkMethods <- function(object){
  errors <- character()
  for(i in seq_along(object)){
    samp_err <- .getMethErr(object[[i]], "colnames")
    if(!is.null(samp_err)){
      errors <- c(errors, paste0("Element [", i, "] of ", samp_err))
    }
    feat_err <- .getMethErr(object[[i]], "rownames")
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

##  Make sure Elist is valid before checking all of the sample names 

S4Vectors::setValidity2("Elist", .validElist)

#' Show method for \code{\linkS4class{Elist}} class
#'
#' @param object A \code{\linkS4class{Elist}} object. 
#' @return Returns a summary of contents for the \code{\linkS4class{Elist}} class
#' exportMethod show
setMethod("show", "Elist", function(object){
		  o_class <- class(object)
		  elem_cl <- vapply(object, class, character(1))
		  o_len <- length(object)
		  o_names <- names(object)
		  sampdim <- vapply(object, FUN = function(obj) { length(colnames(obj)) }, FUN.VALUE = integer(1))
		  featdim <- vapply(object, FUN = function(obj) { length(rownames(obj)) }, FUN.VALUE = integer(1))
		  cat(sprintf('"%s"', o_class), "class object of length", paste0(o_len, ':'),
			  sprintf('\n [%i] %s: "%s" - %s samples, %s features', seq(o_len), o_names, elem_cl, sampdim, featdim), "\n") 
		  })

