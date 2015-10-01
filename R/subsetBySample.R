#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}} object 
#' @param logicID A list of logical indicators to be used for subsetting
subsetBySample <- function(MultiAssay, logicID){
  newPheno <- .subPheno(MultiAssay, .getIndexLogical(MultiAssay, idlist))
  newMap <- Map(subset, MultiAssay@sampleMap, idlist)
  newSubset <- Map(subsetSample, MultiAssay@elist, sapply(newMap, "[", 2))
  # Clone and replace slot? 
  MultiAssay@sampleMap <- newMap
  MultiAssay@elist <- S4Vectors::SimpleList(newSubset)
  MultiAssay@masterPheno <- newPheno
  return(MultiAssay)
}
