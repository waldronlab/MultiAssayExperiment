#' Convert map from data.frame or DataFrame to list
#' 
#' \code{mapToList} returns a \code{list} of DataFrames 
#' 
#' @param mapdf A \code{data.frame} or \code{DataFrame} object with
#' identifiers in the first column
#' @param assayCol A character vector of length one indicating the assay
#' names column 
#' @return A \code{list} object of DataFrames for each assay
#' @export mapToList
mapToList <- function(mapdf, assayCol = NULL) {
    if (is.null(assayCol)) {
        stop("Provide assaynames column reference")
    }
    if (!assayCol %in% colnames(mapdf)) {
        stop("assayname column not found in dataframe")
    }
    if (!is.character(assayCol) | length(assayCol) != 1L) {
        stop("assayname must be a string")
    }
    if (is(mapdf, "data.frame")) {
        mapdf <- S4Vectors::DataFrame(mapdf)
    }
    newList <- S4Vectors::split(mapdf, mapdf[, assayCol])
    newList <- lapply(newList, "[", c(1, 2))
    return(newList)
} 
