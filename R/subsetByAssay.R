#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what assay(s) to select  
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByAssay <- function(MultiAssayExperiment, assayIndicator) {
  newSubset <- Elist(MultiAssayExperiment)[assayIndicator]
  listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
  newMap <- listMap[assayIndicator]
  newMap <- .convertList(newMap)
  sampleMap(MultiAssayExperiment) <- newMap
  Elist(MultiAssayExperiment) <- newSubset
  return(MultiAssayExperiment)
}
