#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a \code{numeric}, \code{character} or \code{logical} or \code{\link{MultiAssayView}} object indicating what assay(s) to select  
#' @param drop logical (default FALSE) whether to return a \code{list} of selected experiments
#' @return A \code{\link{MultiAssayExperiment}} object or \code{list} if drop paramater is set to TRUE
subsetByAssay <- function(MultiAssayExperiment, assayIndicator){
  if(!is(assayIndicator, "MultiAssayView")){
    assayIndicator <- MultiAssayView(MultiAssayExperiment, assayIndicator, "assays")  
  } else {
    if((assayIndicator@type != "assays")){
      stop("Provide an assay type MultiAssayView class")
    }
  }
  assayDrops <- .convertList(assayIndicator@drops, type = "assays")
  assayIndicator <- unlist(assayIndicator@keeps[match(names(MultiAssayExperiment), names(assayIndicator))])
  listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
  newMap <- listMap[assayIndicator]
  newMap <- .convertList(newMap)
  newSubset <- Elist(MultiAssayExperiment)[assayIndicator]
    MultiAssayExperiment@sampleMap <- newMap
    MultiAssayExperiment@Elist <- newSubset
    MultiAssayExperiment@drops <- c(MultiAssayExperiment@drops, list(assayDrops))
  return(MultiAssayExperiment)
}
