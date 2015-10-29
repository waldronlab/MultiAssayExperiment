#' Identify by feature
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param feature A \code{character} vector containing feature names
#' @return A logical list of matched assays
#' @export identifyByFeature
identifyByFeature <- function(MultiAssay, feature, requireall = FALSE){
	featResults <- lapply(MultiAssay@elist, features)
	logicresult <- sapply(seq(length(MultiAssay)),
						  function(i, feats) {feats %in% unlist(featResults[i]) }, feats = feature)
	if(!is.vector(logicresult)){
		logicvector <- apply(logicresult, 2, ifelse(requireall, all, any)) 
	} else {
		logicvector <- logicresult
	}
	logiclist <- sapply(seq(length(featResults)), 
						function(i, feats) {unlist(featResults[i]) %in% feats}, feats = feature)
	logiclist <- lapply(logiclist, `!`)

	dropped <- Map("[", featResults, logiclist)

	newIdentify <- new("Identify", 
					   logreturn = logicvector,
					   drops = dropped)
	return(newIdentify)
}
