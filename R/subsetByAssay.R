#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a \code{numeric}, \code{character} or \code{logical} or \code{\link{MultiAssayView}} object indicating what assay(s) to select  
#' @param drop logical (default FALSE) Indicates whether to return a \code{list} of selected experiments
#' @return A \code{\link{MultiAssayExperiment}} object or \code{list} if drop paramater is set to TRUE
#' @export subsetByAssay
subsetByAssay <- function(MultiAssay, assayIndicator, drop = FALSE){
  if(!is(assayIndicator, "MultiAssayView")){
    assayIndicator <- MultiAssayView(MultiAssay, assayIndicator, "assays")  
  } else {
    if((assayIndicator@type != "assays")){
      stop("Provide an assay type MultiAssayView class")
    }
  }
  assayDrops <- .convertList(assayIndicator@drops, type = "assays")
  assayIndicator <- unlist(assayIndicator@keeps[match(names(MultiAssay), names(assayIndicator))])
  listMap <- toListMap(sampleMap(MultiAssay), "assayname")
  newMap <- listMap[assayIndicator]
  newMap <- .convertList(newMap)
  newSubset <- Elist(MultiAssay)[assayIndicator]
  if(drop){
      return(as.list(newSubset))
  } else {
    MultiAssay@sampleMap <- newMap
    MultiAssay@Elist <- newSubset
    MultiAssay@drops <- c(MultiAssay@drops, list(assayDrops))
  }
  return(MultiAssay)
}
