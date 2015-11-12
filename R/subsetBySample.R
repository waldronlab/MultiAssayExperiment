
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
	newMap <- .ldmap(identifier@keeps)
	chars <- lapply(identifier@keeps, function(x) as.character(unlist(x[,2])))
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
