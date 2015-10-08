#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param logicID A list of logical indicators to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, logicID){
  newPheno <- .subPheno(MultiAssay, .getIndexLogical(MultiAssay, logicID))
  newMap <- Map(subset, MultiAssay@sampleMap, logicID)
  newSubset <- suppressWarnings(S4Vectors::mendoapply(subsetSample, MultiAssay@elist, sapply(newMap, "[", 2)))
  # Clone or replace method for slot??
  MultiAssay@sampleMap <- newMap
  MultiAssay@elist <- newSubset
  MultiAssay@masterPheno <- newPheno
  return(MultiAssay)
}
