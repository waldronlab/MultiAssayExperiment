#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param featureIndicator Either a \code{character} vector or \code{\linkS4class{MultiAssayView}} class object containing feature names
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByRow <- function(MultiAssay, featureIndicator, drop = FALSE, ...){ 
  if(!is(featureIndicator, "MultiAssayView")){
    featureIndicator <- MultiAssayView(MultiAssay, featureIndicator, "rownames")
  } else {
    if(featureIndicator@type != "rownames"){
      stop("Provide a rownames type MultiAssayView class")
    }
  }
  MultiAssay@drops <- c(MultiAssay@drops, list(.convertList(featureIndicator, "drops")))
  featureMaker <- lapply(featureIndicator@keeps, unlist)
  # MultiAssay@Elist <- Elist(MultiAssay)[featureMaker]
  MultiAssay@Elist <- mapply(function(x, y, i, ...){x[y, , drop = i, ...]}, x = Elist(MultiAssay), y = featureMaker, i = drop, ...)
  #  MultiAssay@Elist <- mendoapply(Elist(MultiAssay), subsetFeature, featureMaker)
  return(MultiAssay)
}
