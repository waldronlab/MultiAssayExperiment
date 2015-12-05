#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param identifier A \linkS4class{stage} class object to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, identifier){
	if(identifier@type != "samples"){
		stop("stage class should be of samples!")
	}
	newMap <- getMap(identifier)
	subsetor <- lapply(identifier@keeps, function(x) unlist(x[,2]))
	##  mendoapply not working
	## 	newSubset <- S4Vectors::mendoapply(subsetSample, MultiAssay@elist, chars) 
	newSubset <- mapply(subsetSample, MultiAssay@elist, subsetor)
	newSubset <- elist(newSubset)
	# Clone or replace method for slot??
	MultiAssay@sampleMap <- newMap
	MultiAssay@elist <- newSubset
	MultiAssay@drops <- .convertList(identifier, "drops")
	return(MultiAssay)
}
