#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a \code{numeric}, \code{character} or \code{logical} object indicating what assay(s) to select  
#' @param drop logical (default FALSE) whether to return a \code{list} of selected experiments
#' @return A \code{\link{MultiAssayExperiment}} object or \code{list} if drop paramater is set to TRUE
subsetByAssay <- function(MultiAssayExperiment, assayIndicator){
  newSubset <- Elist(MultiAssayExperiment)[assayIndicator]
  listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
  newMap <- listMap[assayIndicator]
  newMap <- .convertList(newMap)
  sampleMap(MultiAssayExperiment) <- newMap
  Elist(MultiAssayExperiment) <- newSubset
  return(MultiAssayExperiment)
}
