#' Subset MultiAssayExperiment object
#' \code{subset} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}} object
#' @param indicator A \code{logical} or \code{character} vector or \code{Identify} class object to use for subsetting
#' @param method A \code{character} vector of length one designating to subset either by samples, features, or assays
#' @describeIn MultiAssayExperiment
#' export subset
setGeneric("subset", function(x, indicator, ...) standardGeneric("subset"))
setMethod("subset", signature("MultiAssayExperiment", "ANY"), function(x, indicator, ...){
		  if(method == "samples"){
			  indicator <- identifyBySample(MultiAssay = x, indicator)
			  return(subsetBySample(MultiAssay = x, indicator))
		  } else if(method == "features"){
			  indicator <- identifyByFeature(MultiAssay = x, indicator)
			  return(subsetByFeature(MultiAssay = x, indicator))
		  } else if(method == "assays"){
			  return(subsetByAssay(MultiAssay = x, indicator, ...))
		  } else {
			  stop("Please provide a valid subsetting 'by' argument!")
		  }
})

setMethod("subset", signature("MultiAssayExperiment", "Identify"), function(x, indicator, ...){
		  return(subset(MultiAssay = x, indicator))
})
