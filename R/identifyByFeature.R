#' Identify by feature
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param feature A \code{character} vector containing feature names
#' @return A logical list of matched assays
#' @export identifyByFeature
identifyByFeature <- function(MultiAssay, feature, requireall = FALSE){
	featResults <- lapply(MultiAssay@elist, featExtractor)
	logicMatrix <- sapply(seq(length(MultiAssay)),
						  function(i, feats) {feats %in% unlist(featResults[i]) }, feats = feature)
	if(!is.vector(logicMatrix)){
		if(requireall){
			return(apply(logicMatrix, 2, all))
		}else{
			return(apply(logicMatrix, 2, any))
		}
	}else{
		return(logicMatrix)
	}
}
