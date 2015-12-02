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
  v_len <- vapply(object@keeps, FUN = function(x) {nrow(x)}, FUN.VALUE = integer(1))
  cat("A", sprintf('"%s"', o_class), "class object of length", paste0(o_len, ':'),
      "\nIdentifiers: ")
  cat(o_ids, sep = ", ")
  cat("\n Staged by: ", '"', stage_type, '"', sep = "")
  cat(sprintf('\n [%i] %s: %s %s', seq(o_len), o_names, v_len, stage_type), "\n")
})
