.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.getLogicalNames <- function(object, ids){
return(lapply(object@sampleMap, function(map) { map[, 1] %in% ids}))
}

.getNamesLogical <- function(object, logiID){
subList <- Map(subset, object@sampleMap, logiID)
usedNames <- Reduce(union, sapply(subList, "[", 1))
return(rownames(object@masterPheno)[match(usedNames, rownames(object@masterPheno))])
}

.getIndexLogical <- function(object, logiID){
  subList <- Map(subset, object@sampleMap, logiID)
  usedNames <- Reduce(union, sapply(subList, "[", 1))
  return(match(usedNames, rownames(object@masterPheno)))
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
	return(.getLogicalNames(MAobject, iders))
}
