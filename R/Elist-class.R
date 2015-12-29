.getMethErr <- function(object, my_fun){
	obj_cl <- class(object)
	e_class <- class(try(get(my_fun)(object), silent = TRUE))
	if(e_class == "try-error"){
		msg <-  paste0("class ", obj_cl, " should have a '", my_fun, "' method")
		return(msg)
	}
	NULL
}

## Check class conforms to API
.hasMethods <- function(object, my_fun){
  obj_cl <- class(object)
  if(any(my_fun %in% c("[", "assay"))){
    if(any(obj_cl %in% c("GRangesList", "RangedRaggedAssay"))){
      return(hasMethod(my_fun, signature = c(obj_cl, "ANY"),
                where = c("package:GenomicRanges", "package:IRanges", "package:SummarizedExperiment")))
    } else if(is(object, "RangedSummarizedExperiment")){
      return(hasMethod(my_fun, signature = c(class(object), "missing"), 
                where = c("package:GenomicRanges", "package:IRanges", "package:SummarizedExperiment")))
    } else {
      return(hasMethod(my_fun, signature = c(obj_cl, "ANY")))
    }
  }
  return(hasMethod(my_fun, signature = obj_cl))
}

.getNameErr <- function(object){
  if(inherits(object, "RangedRaggedAssay")){
    if(is.null(names(object))){
      msg <- paste("names in", obj_cl, "are NULL")
      return(msg)
    } else { NULL } 
  } else { NULL }
}

### ==============================================
### Elist class
### ----------------------------------------------

#' A container for multi-experiment data
#' 
#' The \code{Elist} class is a container that builds on 
#' the \code{\link[S4Vectors]{SimpleList-class}} with additional 
#' checks for consistency in experiment names and length. It contains a \code{SimpleList}
#' of experiments with sample identifiers. One element present per experiment performed.  
#' 
#' @inheritParams S4Vectors::SimpleList
#' @exportClass Elist
.Elist <- setClass("Elist", contains = "SimpleList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

setGeneric("Elist", function(x) standardGeneric("Elist"))

#' Convert from \code{\link[S4Vectors]{SimpleList-class}}
#' to the multi-experiment data container
#'
#' @param x A \code{list} object
#' @return An \code{\linkS4class{Elist}} class object
#' @exportMethod Elist
setMethod("Elist", "ANY", function(x){
  .Elist(S4Vectors::SimpleList(x))
})
setMethod("Elist", "missing", function(x){
 .Elist(S4Vectors::SimpleList(list())) 
})

##
## Validity ---------------------------------
##

.checkMethodsTable <- function(object){
  supportedMethods <- c("colnames", "rownames", "[", "assay")
  errors <- character(length(object)*length(supportedMethods))
  for(i in seq_along(object)){
    for(j in seq_along(supportedMethods)){
      if(!(.hasMethods(object[[i]], supportedMethods[[j]]))){
        errors <- c(errors, paste0("Element [", i, "] of class ", class(object[[i]]), " does not have a ", supportedMethods[[j]], " method"))
      }
    }
  }
  return(errors)
}

.forceEvalMethods <- function(object){
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
    msg <- "Non-unique names provided"
    errors <- c(errors, msg)
  } 
  if(length(errors) == 0L){
    NULL
  } else { errors }
}

.validElist <- function(object){
  if(length(object) != 0L){
  c(.forceEvalMethods(object),
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

