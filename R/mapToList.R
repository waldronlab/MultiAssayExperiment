#' Convert map from data.frame or DataFrame to list and vice versa
#' 
#' The \code{mapToList} function provides a convenient way of reordering a
#' \code{data.frame} to a \code{list}. The \code{listToMap} function does the 
#' opposite by taking a \code{list} and converting it to \code{DataFrame}.
#' 
#' @param dfmap A \code{data.frame} or \code{DataFrame} object with
#' identifiers in the first column
#' @param assayCol A character vector of length one indicating the assay
#' names column 
#' @return A \code{list} object of DataFrames for each assay
#' @aliases listToMap 
#' @example inst/scripts/listToMap-Ex.R
#' @export mapToList
mapToList <- function(dfmap, assayCol = "assay") {
    if (!S4Vectors::isSingleString(assayCol))
        stop("assay column name must be a single string")
    if (!assayCol %in% colnames(dfmap))
        stop("assay column not found in dataframe")
    assayColIndex <- which(names(dfmap) == assayCol)
    if (inherits(dfmap, "data.frame")) {
        dfmap <- S4Vectors::DataFrame(dfmap)
    }
    assayOrder <- unique(dfmap[[assayCol]])
    newList <- S4Vectors::split(dfmap[, -assayColIndex], dfmap[[assayCol]])
    ## Preserve the order of the assaynames in the map!
    newList <- newList[assayOrder]
    return(newList)
} 
