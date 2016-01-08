#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param rowIndicator A \code{character} vector class object containing feature names
#' @param drop logical (default FALSE) whether to coerce lowest possible dimension after subsetting
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssayExperiment, rowIndicator, ...) { 
  if(is(rowIndicator, "character")){
  loglistmatch <- lapply(rownames(MultiAssayExperiment), function(charElem) {
    charElem %in% rowIndicator })
  } else {
    stop("rowIndicator must be character")
  }
  # MultiAssayExperiment@Elist <- Elist(MultiAssayExperiment)[featureMaker]
  MultiAssayExperiment@Elist <- Elist(mapply(function(x, i, j, ..., drop = FALSE) 
    {x[i, , ..., drop = FALSE]}, x = Elist(MultiAssayExperiment), i = loglistmatch, ...))
  #  MultiAssayExperiment@Elist <- mendoapply(Elist(MultiAssayExperiment), subsetFeature, featureMaker)
  return(MultiAssayExperiment)
}
