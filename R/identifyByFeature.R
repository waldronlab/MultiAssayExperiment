#' Identify by feature
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param feature A \code{character} vector containing feature names
#' @return A logical list of matched assays
#' @export identifyByFeature
identifyByFeature <- function(MultiAssay, feature){
	featResults <- lapply(MultiAssay@elist, featExtractor)
	return(sapply(seq(length(MultiAssay)),
				  function(i, feats) {feats %in% unlist(featResults[i]) }, feats = feature))
}

