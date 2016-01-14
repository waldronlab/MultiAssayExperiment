#' Subset MultiAssayExperiment object
#' \code{subsetByColumn} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object 
#' @param colIndicator A \linkS4class{MultiAssayView} class object to be used for subsetting
#' @param drop logical (default FALSE) whether to coerce lowest possible dimension after subsetting
subsetByColumn <- function(MultiAssayExperiment, colIndicator){
 if(is.character(colIndicator)){
  loglistmatch <- lapply(colnames(MultiAssayExperiment), function(charElem) {
    charElem %in% colIndicator})
  if(!any(unlist(loglistmatch))){
    stop("No matches found")
  }
 } else {
   stop("colIndicator must be character")
 }
  ## 	newSubset <- S4Vectors::mendoapply(subsetSample, Elist(MultiAssayExperiment), subsetor) 
  newSubset <- mapply(function(x, i, j, drop) {x[, j, drop = FALSE]}, x = Elist(MultiAssayExperiment), j = loglistmatch)
  newSubset <- Elist(newSubset)
  listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
  logmapInd <- lapply(listMap, function(x) { x[,2] %in% colIndicator })
  newMap <- mapply(function(x, y) {x[y,]}, listMap, logmapInd)
  newMap <- Filter(function(x) {!isEmpty(x)}, newMap)
  newMap <- .convertList(newMap)
  # Clone or replace method for slot??
  sampleMap(MultiAssayExperiment) <- newMap
  Elist(MultiAssayExperiment) <- newSubset
# MultiAssayExperiment@drops <- c(MultiAssayExperiment@drops, list(.convertList(colIndicator, "drops")))
  return(MultiAssayExperiment)
}
