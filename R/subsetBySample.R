#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param logic_cl An \linkS4class{Identify} class object to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, logic_cl){
	if(logic_cl@type != "samples"){
		stop("Identify class should be of samples!")
	}
	phenoIndex <- match(logic_cl@query, rownames(MultiAssay@masterPheno))
	newPheno <- .subPheno(MultiAssay, phenoIndex)
	listMap <- toListMap(MultiAssay@sampleMap, "assayname")
	listMap <- listMap[order(names(logic_cl@keeps))]
	newMap <- Map(subset, listMap, logic_cl@keeps)
	chars <- lapply(newMap, function(x) as.character(unlist(x[2])))
	newSubset <- S4Vectors::mendoapply(subsetSample, MultiAssay@elist, chars) 
	# Clone or replace method for slot??
	MultiAssay@sampleMap <- newMap
	MultiAssay@elist <- newSubset
	MultiAssay@masterPheno <- newPheno
	return(MultiAssay)
}
