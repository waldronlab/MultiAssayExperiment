#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param featureIndicator Either a \code{character} vector or \code{\linkS4class{MultiAssayView}} class object containing feature names
#' @param drop logical indicating whether to coerce to lowest possible dimension after subsetting
#' @param ... Additional arguments to pass to low level subsetting function
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssayExperiment, featureIndicator, drop = FALSE, ...){ 
  if(!is(featureIndicator, "MultiAssayView")){
    featureIndicator <- MultiAssayView(MultiAssayExperiment, featureIndicator, "rownames")
  } else {
    if(featureIndicator@type != "rownames"){
      stop("Provide a rownames type MultiAssayView class")
    }
  }
  MultiAssayExperiment@drops <- c(MultiAssayExperiment@drops, list(.convertList(featureIndicator, "drops")))
  featureMaker <- lapply(featureIndicator@keeps, unlist)
  # MultiAssayExperiment@Elist <- Elist(MultiAssayExperiment)[featureMaker]
  MultiAssayExperiment@Elist <- mapply(function(x, y, i, ...){x[y, , drop = i, ...]}, x = Elist(MultiAssayExperiment), y = featureMaker, i = drop, ...)
  #  MultiAssayExperiment@Elist <- mendoapply(Elist(MultiAssayExperiment), subsetFeature, featureMaker)
  return(MultiAssayExperiment)
}
