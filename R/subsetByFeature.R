#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param featureIndicator Either a \code{character} vector or \code{\linkS4class{Stage}} class object containing feature names
#' @return A \code{\link{MultiAssayExperiment}} object 
#' @export subsetByFeature
subsetByFeature <- function(MultiAssay, featureIndicator, ...){ 
  if(!is(featureIndicator, "Stage")){
    featureIndicator <- Stage(MultiAssay, featureIndicator, "features")
  } else {
    if(featureIndicator@type != "features"){
    stop("Provide a feature type Stage class!")
  }
  }
  if(length(MultiAssay@drops)==0L){
    sequence <- paste0("Operation_", 1)
    MultiAssay@drops <- list(.convertList(featureIndicator, "drops"))
    names(MultiAssay@drops) <- sequence
  } else {
    num <- as.numeric(gsub("Operation_", "", names(MultiAssay@drops)[length(MultiAssay@drops)]))+1
    sequence <- paste0("Operation_", num)
    newElement <- list(.convertList(identifier, "drops"))
    names(newElement) <- sequence
    MultiAssay@drops <- c(MultiAssay@drops, newElement)
  }
  featureMaker <- lapply(featureIndicator@keeps, unlist)
  MultiAssay@Elist <- Elist(Map(subsetFeature, MultiAssay@Elist, featureMaker))
#  MultiAssay@Elist <- mendoapply(MultiAssay@Elist, subsetFeature, featureMaker)
	return(MultiAssay)
}

