.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.getLogicalNames <- function(object, ids){
logn <- lapply(object@sampleMap, function(map) { 
   completemap <- map[!is.na(map[, 2]),]
   completemap[, 1] %in% ids })
return(logn)
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
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}}
#' @param j A numeric vector referencing masterPheno data.frame row numbers
#' @return A logical list of matched sample references
#' @export identifyBySample
identifyBySample <- function(MultiAssay, j){
	iders <- rownames(.subPheno(MultiAssay, j))
	return(.getLogicalNames(MultiAssay, iders))
}
