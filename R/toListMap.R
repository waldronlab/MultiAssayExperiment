#' Convert data.frame map to list map
#' \code{toListMap} returns a \code{list} 
#' 
#' @param dfmap A data frame object with identifiers in the first column
#' @param assayCol A character vector of length one indicating the assay names column 
#' @export
toListMap <- function(dfmap, assayCol = NULL){
	if(is.null(assayCol)){
		stop("Provide the name of the column for assay names!")
	}
	if(!assayCol %in% colnames(dfmap)){
		stop("Assay names column not found in dataframe!")
	}
	if(!is.character(assayCol) | length(assayCol) != 1L){
		stop("Assay name must be a string!")
	}
	colIndex <- which(colnames(dfmap) == assayCol)
	newList <- split(dfmap[-colIndex], dfmap[assayCol])
	return(lapply(newList, FUN = function(x) { rownames(x) <- NULL
				  x} ))
}
