#' Convert data.frame map to list map
#' \code{toListMap} returns a \code{list} 
#' 
#' @param dfmap A data frame object with identifiers in the first column
#' @export
toListMap <- function(dfmap){
	newlist <- lapply(seq_along(dfmap)[-1], FUN = function(col, mydf) { 
					  cbind(mydf[1], mydf[col]) }, mydf = dfmap)
	names(newlist) <- names(dfmap)[-1]
	return(newlist)
}
