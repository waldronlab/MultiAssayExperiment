#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param identifier A \linkS4class{stage} class object to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, identifier){
	if(is(identifier, "stage") && getElement(identifier, "type") != "samples"){
		stop("stage class should be of samples!")
	} else {
	newMap <- getMap(identifier)
	subsetor <- lapply(identifier@keeps, function(x) unlist(x[,2]))
	}
	##  mendoapply not working
	## 	newSubset <- S4Vectors::mendoapply(subsetSample, MultiAssay@elist, chars) 
	newSubset <- mapply(subsetSample, MultiAssay@elist, subsetor)
	newSubset <- elist(newSubset)
	# Clone or replace method for slot??
	MultiAssay@sampleMap <- newMap
	MultiAssay@elist <- newSubset
	if(length(MultiAssay@drops)==0L){
	  sequence <- paste0("Operation_", 1)
	  MultiAssay@drops <- list(.convertList(identifier, "drops"))
	  names(MultiAssay@drops) <- sequence
	} else {
	  num <- as.numeric(gsub("Operation_", "", names(MultiAssay@drops)[length(MultiAssay@drops)]))+1
	  sequence <- paste0("Operation_", num)
	  newElement <- list(.convertList(identifier, "drops"))
	  names(newElement) <- sequence
	  MultiAssay@drops <- c(MultiAssay@drops, newElement)
	}
	return(MultiAssay)
}
