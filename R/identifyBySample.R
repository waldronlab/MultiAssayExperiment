.arrangeMap <- function(map, ordering){
	newOrd <- do.call(c, lapply(ordering, function(ord) { 
					 which(map[, 1] == ord) } ))
	return(newOrd)
}

.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.cutMap <- function(map, lid){
	compdf <- map[!is.na(map[, 2]), ]
	return(compdf)
}

.separateMap <- function(object, ids){
	DFsampleMap <- S4Vectors::DataFrame(object@sampleMap)
	listDFsampleMap <- toListMap(DFsampleMap, "assayname")
	listDFsampleMap <- listDFsampleMap[order(names(object@elist))]
	loglistmatch <- lapply(listDFsampleMap, function(map) { map[,"master"] %in% ids })
	keeps <- Map(function(x, y) { x[y,] }, listDFsampleMap, loglistmatch)
	orderIndex <- lapply(keeps, .arrangeMap, ids)
	orderedKeeps <- Map(function(x, y) { x[y, ] }, keeps, orderIndex)
	duoMap <- list(keeps = orderedKeeps,
				   drops = Map(function(x, y) { x[y,] }, listDFsampleMap, lapply(loglistmatch, `!`)))
	return(duoMap)
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
	newIdentify <- new("Identify",
					   query = iders,
					   keeps = biMap[["keeps"]],
					   drops = biMap[["drops"]],
					   type = "samples")
	return(newIdentify)
}
