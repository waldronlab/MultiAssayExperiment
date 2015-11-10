.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.cutMap <- function(map, lid){
	compdf <- map[!is.na(map[, 2]), ]
	return(compdf)
}

.getLogicalNames <- function(object, ids){
	listMap <- toListMap(object@sampleMap, "assayname")
	logn <- lapply(listMap, function(map) { map[, 1] %in% ids } )
	return(logn)
}

.getNamesLogical <- function(object, logiID){
	listMap <- toListMap(object@sampleMap, "assayname")
	subList <- Map(subset, listMap, logiID)
	usedNames <- Reduce(union, sapply(subList, "[", 1))
	return(rownames(object@masterPheno)[match(usedNames, rownames(object@masterPheno))])
}

.getIndexLogical <- function(object, logiID){
	listMap <- toListMap(object@sampleMap, "assayname")
	listMap <- listMap[order(names(logiID@keeps))]
	subList <- Map(subset, listMap, logiID@keeps)
	usedNames <- logiID@query
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
	sampResults <- samples(MultiAssay)
	iders <- rownames(.subPheno(MultiAssay, j))
	logiclist <- .getLogicalNames(MultiAssay, iders)
	logiclist <- logiclist[order(names(sampResults))]
	revlogResult <- lapply(logiclist, "!")
	dropped <- Map("[", sampResults, revlogResult)
	newIdentify <- new("Identify",
					   query = iders,
					   keeps = logiclist,
					   drops = dropped,
					   type = "samples")
	return(newIdentify)
}
