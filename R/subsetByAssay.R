#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets,
#' fuzzy matching is possible but exact names recommended
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param assayName Either a numeric or logical vector indicating what assay(s) to select  
#' @param drop logical Indicates whether to drop the unmatched assays from the object
#' @return A \code{\link{MultiAssayExperiment}} object or other if drop paramater is set to TRUE
#' @export subsetByAssay
subsetByAssay <- function(MultiAssay, assayName, drop = FALSE){
	lengthNames <- c(seq(length(MultiAssay)), names(MultiAssay))
if(!is.numeric(assayName)){
assayName <- vapply(as.character(assayName),
	FUN = function(x){ grep(x, lengthNames, ignore.case = TRUE, value = TRUE) },
	FUN.VALUE = character(1))
}
	if(all(assayName %in% lengthNames)){
		newMap <- MultiAssay@sampleMap[assayName]
		newSubset <- MultiAssay@elist[assayName]
		if(length(assayName)==1 & drop){
			return(as.list(newSubset))
		} else {
			MultiAssay@sampleMap <- newMap
			MultiAssay@elist <- newSubset
			return(MultiAssay)
		}
	} else {
		stop("Given assay vector is not in the assay list!")
	}
}
