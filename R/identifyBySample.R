.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.getSubsetLogical <- function(object, ids){
return(lapply(object@sampleMap, function(map) { map[, 1] %in% ids}))
}

#' Identify samples corresponding to row index in masterPheno
#' \code{identifyBySample} returns a \code{logical} 
#' 
#' This function uses row indices of the masterPheno to return a logical list
#' 
#' @param MAobject A \code{\linkS4class{MultiAssayExperiment}}
#' @param j A numeric vector referencing masterPheno data.frame row numbers
identifyBySample <- function(MAobject, j){
	iders <- rownames(.subPheno(MAobject, j))
	return(.getSubsetLogical(MAobject, iders))
}
