#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param rowIndicator A \code{character} vector or \code{GRanges} class object
#' containing feature names or ranges
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssayExperiment, rowIndicator, ...) {
  if (!is.character(rowIndicator) && !is(rowIndicator, "GRanges")) {
    stop("rowIndicator must be character or GRanges")
  }
  hitList <- getHits(MultiAssayExperiment, rowIndicator, ...)
  Elist(MultiAssayExperiment) <-
    Elist(mapply(function(x, i, j, ..., drop = FALSE) {
      x[i, , ..., drop = FALSE]
    }, x = Elist(MultiAssayExperiment), i = hitList, ...))
  return(MultiAssayExperiment)
} 
