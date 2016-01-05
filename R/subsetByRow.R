#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param rowIndicator Either a \code{character} vector or \code{\linkS4class{MultiAssayView}} class object containing feature names
#' @param drop logical (default FALSE) whether to coerce lowest possible dimension after subsetting
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssayExperiment, rowIndicator, drop = FALSE, ...){ 
  if(!is(rowIndicator, "MultiAssayView")){
    rowIndicator <- MultiAssayView(MultiAssayExperiment, rowIndicator, "rownames")
  } else {
    if(rowIndicator@type != "rownames"){
      stop("Provide a rownames type MultiAssayView class")
    }
  }
  MultiAssayExperiment@drops <- c(MultiAssayExperiment@drops, list(.convertList(rowIndicator, "drops")))
  featureMaker <- lapply(rowIndicator@keeps, unlist)
  # MultiAssayExperiment@Elist <- Elist(MultiAssayExperiment)[featureMaker]
  MultiAssayExperiment@Elist <- Elist(mapply(function(x, i, j, ..., drop){x[i, , ..., drop = drop]}, x = Elist(MultiAssayExperiment), i = featureMaker, drop = drop, ...))
  #  MultiAssayExperiment@Elist <- mendoapply(Elist(MultiAssayExperiment), subsetFeature, featureMaker)
  return(MultiAssayExperiment)
}
