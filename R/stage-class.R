.getVLen <- function(obj){
  list_ob <- S4Vectors::split(obj@keeps, obj@keeps[["assayname"]])
  browser()
  sapply(seq_along(list_ob), function(i, x) {
    if(is.na(x[[i]][["feature"]])){
      0L
    } else {
      nrow(x[[i]]) 
    }
  }, x = list_ob)
}

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

#' Show method for \code{\linkS4class{stage}} class
#' 
#' @param object A \code{\linkS4class{stage}} class object
#' @return A summary of \code{\linkS4class{stage}} class contents 
#' @exportMethod show
setMethod("show", "stage", function(object){
  o_class <- class(object)
  o_len <- length(object)
  o_names <- names(object)
  stage_type <- type(object)
  o_ids <- query(object)
  #  v_len <- .getVLen(object)
  if(stage_type != "assays"){
    my_fun <- nrow
  } else {
    my_fun <- function(logic){
      if(logic) "keep"
      else "drop"
    }
  }
  v_ops <- sapply(object@keeps, FUN = function(x) {my_fun(x)})
  cat("A", sprintf('"%s"', o_class), "class object of length", paste0(o_len, ':'),
      "\nQuery: ")
  cat(o_ids, sep = ", ")
  cat("\n Staged by: ", '"', stage_type, '"', sep = "")
  cat(sprintf('\n [%i] %s: %s', seq(o_len), o_names, paste(v_ops, if(is.numeric(v_ops)){ifelse(v_ops == 1L, gsub("s$", "", stage_type), stage_type)}), "\n"))
})
