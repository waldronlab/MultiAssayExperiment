### ==============================================
### Stage class 
### ==============================================

#' An identifier class used for staging a subset operation
#' 
#' @slot query Any class indicator needed to subset
#' @slot keeps A \code{list} indicating valid matches in each assay
#' @slot drops A \code{list} of excluded information due to subsetting
#' @slot type A \code{character} vector indicating method used to search
#' @exportClass Stage
setClass("Stage", 
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

S4Vectors::setValidity2("Stage", .validStage)

#' Show method for \code{\linkS4class{Stage}} class
#' 
#' @param object A \code{\linkS4class{Stage}} class object
#' @return A summary of \code{\linkS4class{Stage}} class contents 
#' @exportMethod show
setMethod("show", "Stage", function(object){
  o_class <- class(object)
  o_len <- length(object)
  o_names <- names(object)
  stage_type <- type(object)
  o_ids <- rownames(query(object))
  if(stage_type != "assays"){
    my_fun <- function(x) length(na.omit(x[, 1]))
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
