.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.cutMap <- function(map, lid){
	compdf <- map[!is.na(map[, 2]), ]
	return(compdf)
}

.getLogicalNames <- function(object, ids){
	trimMap <- lapply(object@sampleMap, .cutMap)
	logn <- lapply(trimMap, function(map) { map[, 1] %in% ids } )
	return(logn)
}

.getNamesLogical <- function(object, logiID){
	trimMap <- lapply(object@sampleMap, .cutMap)
	subList <- Map(subset, trimMap, logiID)
	usedNames <- Reduce(union, sapply(subList, "[", 1))
	return(rownames(object@masterPheno)[match(usedNames, rownames(object@masterPheno))])
}

.getIndexLogical <- function(object, logiID){
	trimMap <- lapply(object@sampleMap, .cutMap)
	subList <- Map(subset, trimMap, logiID)
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
	sampResults <- lapply(MultiAssay@elist, samples)
	iders <- rownames(.subPheno(MultiAssay, j))
	logiclist <- .getLogicalNames(MultiAssay, iders)
	revlogResult <- lapply(logiclist, "!")
	dropped <- Map("[", sampResults, revlogResult)
	newIdentify <- new("Identify",
					   indim = logiclist, 
					   identifier = iders, 
					   drops = dropped)
	return(newIdentify)
}
