#' Subset MultiAssayExperiment object
#' \code{subset} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}} object
#' @param indicator A \code{logical} vector or object of class \code{Identify} to use for subsetting
#' @param by A \code{character} vector of length one designating to subset either by samples, features, or assays
setMethod("subset", "MultiAssayExperiment", function(x, indicator, by, ...){
		  if(by == "samples"){
			  idclass <- identifyBySample(x, indicator)
			  return(subsetBySample(x, idclass))
		  } else if(by == "features"){
			  idclass <- identifyByFeatures(x, indicator)
			  return(subsetByFeatures(x, idclass))
		  }
})
