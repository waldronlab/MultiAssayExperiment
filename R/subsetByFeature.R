#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param featureIndicator Either a \code{character} vector or \code{\linkS4class{Stage}} class object containing feature names
#' @return A \code{\link{MultiAssayExperiment}} object 
#' @export subsetByFeature
subsetByFeature <- function(MultiAssay, featureIndicator, ...){ 
  if(!is(featureIndicator, "MultiAssayView")){
    featureIndicator <- Stage(MultiAssay, featureIndicator, "rownames")
  } else {
    if(featureIndicator@type != "rownames"){
    stop("Provide a rownames type MultiAssayView class!")
  }
  }
  MultiAssay@drops <- c(MultiAssay@drops, list(.convertList(featureIndicator, "drops")))
  featureMaker <- lapply(featureIndicator@keeps, unlist)
  MultiAssay@Elist <- MultiAssay@Elist[featureMaker]
  #  MultiAssay@Elist <- mendoapply(MultiAssay@Elist, subsetFeature, featureMaker)
	return(MultiAssay)
}

