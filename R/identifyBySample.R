.arrangeMap <- function(map, ordering){
	newOrd <- do.call(c, lapply(ordering, function(ord) { 
					 which(map[, 1] == ord) } ))
	return(newOrd)
}

.getLogicalNames <- function(object, ids){
	listMap <- toListMap(object@sampleMap, "assayname")
	logn <- lapply(listMap, function(map) { map[, 1] %in% ids } )
	return(logn)
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
	if(is.numeric(j)){
		iders <- rownames(.subPheno(MultiAssay, j))
	} else {
		iders <- .subPheno(MultiAssay, j)
	}
	biMap <- .separateMap(MultiAssay, iders)
	newStage <- new("stage",
					   query = iders,
					   keeps = biMap[["keeps"]],
					   drops = biMap[["drops"]],
					   type = "samples")
	return(newStage)
}
