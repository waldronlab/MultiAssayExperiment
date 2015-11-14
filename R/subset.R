#' Subset MultiAssayExperiment object
#' \code{subset} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}} object
#' @param indicator A \code{logical} or \code{character} vector or \code{stage} class object to use for subsetting
#' @param method A \code{character} vector of length one designating to subset either by samples, features, or assays
#' @describeIn MultiAssayExperiment
#' export subset
setMethod("subset", "MultiAssayExperiment", function(x, indicator, ...){
		  if(length(list(...)) > 0L){
			  stop("invalid subsetting")}
		  method <- indicator@type
		  if(method == "samples"){
			  indicator <- identifyBySample(MultiAssay = x, indicator)
			  return(subsetBySample(MultiAssay = x, indicator))
		  } else if(method == "features"){
			  indicator <- identifyByFeature(MultiAssay = x, indicator)
			  return(subsetByFeature(MultiAssay = x, indicator))
		  } else if(method == "assays"){
			  return(subsetByAssay(MultiAssay = x, indicator))
		  } 
})
