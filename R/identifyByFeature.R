#' Identify by feature
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param feature A \code{character} vector containing feature names
#' @param requireall logical (default FALSE) force matching of all features provided
#' @return A logical list of matched assays
#' @export identifyByFeature
identifyByFeature <- function(MultiAssay, feature, requireall = FALSE){
	featResults <- features(MultiAssay) 
	logicresult <- sapply(seq_along(MultiAssay),
						  function(i, feats) {feats %in% unlist(featResults[i]) }, feats = feature)
	if(!is.vector(logicresult)){
		logicvector <- apply(logicresult, 2, ifelse(requireall, all, any)) 
	} else {
		logicvector <- logicresult
	}
	logiclist <- sapply(seq_along(featResults), 
						function(i, feats) {unlist(featResults[i]) %in% feats}, feats = feature)
	names(logiclist) <- names(featResults)
	revLogicList <- lapply(logiclist, `!`)
	dropped <- Map("[", featResults, revLogicList)
	newIdentify <- new("Identify", 
					   query = feature,
					   keeps = logiclist,
					   drops = dropped,
					   type = "features")
	return(newIdentify)
}
