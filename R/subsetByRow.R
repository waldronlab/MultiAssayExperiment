#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param rowIndicator A \code{character} vector class object containing feature names
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssayExperiment, rowIndicator, ...) { 
  if(is.character(rowIndicator)){
  orderedMatches <- lapply(rownames(MultiAssayExperiment), function(charElem) {
    na.omit(charElem[charElem %in% rowIndicator][order(rowIndicator)]) })
  } else {
    stop("rowIndicator must be character or GRanges")
  }
  # MultiAssayExperiment@Elist <- Elist(MultiAssayExperiment)[featureMaker]
  MultiAssayExperiment@Elist <- Elist(mapply(function(x, i, j, ..., drop = FALSE) 
    {x[i, , ..., drop = FALSE]}, x = Elist(MultiAssayExperiment), i = orderedMatches, ...))
  #  MultiAssayExperiment@Elist <- mendoapply(Elist(MultiAssayExperiment), subsetFeature, featureMaker)
  return(MultiAssayExperiment)
}