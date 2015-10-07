#' Subset MultiAssayExperiment object by Feature
#' 
#' Subset MultiAssayExperiment class by provided feature names
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param feature A \code{character} vector containing feature names
#' @return A \code{\link{MultiAssayExperiment}} object 
#' @export subsetByFeature
subsetByFeature <- function(MultiAssay, feature){
	subsetLogic <- identifyByFeature(MultiAssay, feature)  
	MultiAssay <- subsetByAssay(MultiAssay, subsetLogic, drop = FALSE)
	MultiAssay@elist <- endoapply(MultiAssay@elist, subsetFeature, feature)
	return(MultiAssay)
}
