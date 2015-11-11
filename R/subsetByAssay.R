#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a numeric, character or logical vector indicating what assay(s) to select  
#' @param drop logical (default FALSE) Indicates whether to return a \code{list} of selected experiments
#' @return A \code{\link{MultiAssayExperiment}} object or \code{list} if drop paramater is set to TRUE
#' @export subsetByAssay
subsetByAssay <- function(MultiAssay, assayIndicator, drop = FALSE){
	if(is.logical(assayIndicator)){
		if(length(assayIndicator) != length(MultiAssay@elist)){
			stop("Provide a valid logical vector of equal length!")
		}
	} else if(is.character(assayIndicator)){
		if(!all(assayIndicator %in% names(MultiAssay@elist))){
			stop("Provide a vector of valid assay names!")
		} 
	}
	listMap <- toListMap(MultiAssay@sampleMap, "assayname")
	newMap <- listMap[assayIndicator]
	newMap <- .ldmap(newMap)
	newSubset <- MultiAssay@elist[assayIndicator]
	if(drop){
		return(as.list(newSubset))
	} else {
		MultiAssay@sampleMap <- newMap
		MultiAssay@elist <- newSubset
		return(MultiAssay)
	}
}
