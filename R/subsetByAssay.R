#' Subset MultiAssayExperiment object by Assay type
#' 
#' Select which assay(s) to obtain from available datasets
#' 
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object
#' @param assayIndicator Either a \code{numeric}, \code{character} or \code{logical} or \code{\link{Stage}} object indicating what assay(s) to select  
#' @param drop logical (default FALSE) Indicates whether to return a \code{list} of selected experiments
#' @return A \code{\link{MultiAssayExperiment}} object or \code{list} if drop paramater is set to TRUE
#' @export subsetByAssay
subsetByAssay <- function(MultiAssay, assayIndicator, drop = FALSE){
  if(!is(assayIndicator, "MultiAssayView")){
    assayIndicator <- Stage(MultiAssay, assayIndicator, "assays")  
  } else {
    if((assayIndicator@type != "assays")){
      stop("Provide an assay type Stage class!")
    }
  }
  assayDrops <- .convertList(assayIndicator@drops, type = "assays")
  assayIndicator <- unlist(assayIndicator@keeps[match(names(MultiAssay), names(assayIndicator))])
  listMap <- toListMap(MultiAssay@sampleMap, "assayname")
  newMap <- listMap[assayIndicator]
  newMap <- .convertList(newMap)
  newSubset <- MultiAssay@Elist[assayIndicator]
  if(drop){
    return(as.list(newSubset))
  } else {
    MultiAssay@sampleMap <- newMap
    MultiAssay@Elist <- newSubset
    if(length(MultiAssay@drops)==0L){
      sequence <- paste0("Operation_", 1)
      MultiAssay@drops <- list(assayDrops)
      names(MultiAssay@drops) <- sequence
    } else {
      num <- as.numeric(gsub("Operation_", "", names(MultiAssay@drops)[length(MultiAssay@drops)]))+1
      sequence <- paste0("Operation_", num)
      newElement <- list(assayDrops)
      names(newElement) <- sequence
      MultiAssay@drops <- c(MultiAssay@drops, newElement)
    }
    return(MultiAssay)
  }
}
