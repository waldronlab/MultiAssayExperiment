.arrangeMap <- function(map, ordering){
	newOrd <- do.call(c, lapply(ordering, function(ord) { 
					 which(map[, 1] == ord) } ))
	return(newOrd)
}

#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param identifier An \linkS4class{Identify} class object to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, identifier){
	if(identifier@type != "samples"){
		stop("Identify class should be of samples!")
	}
	newMap <- toListMap(MultiAssay@sampleMap, "assayname")
	newMap <- newMap[order(names(identifier@keeps))]
	newMap <- Map(base::subset, newMap, identifier@keeps)
	orderIndex <- lapply(newMap, .arrangeMap, identifier@query)
	orderedMap <- Map(function(x, y) { x[y, ] }, newMap, orderIndex)
	chars <- lapply(orderedMap, function(x) as.character(unlist(x[2])))
	newMap <- .ldmap(orderedMap)
	##  mendoapply not working
	## 	newSubset <- S4Vectors::mendoapply(subsetSample, MultiAssay@elist, chars) 
	newSubset <- mapply(subsetSample, MultiAssay@elist, chars)
	newSubset <- elist(newSubset)
	# Clone or replace method for slot??
	MultiAssay@sampleMap <- newMap
	MultiAssay@elist <- newSubset
	MultiAssay@drops <- identifier@drops
	return(MultiAssay)
}
