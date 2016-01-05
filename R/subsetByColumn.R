#' Subset MultiAssayExperiment object
#' \code{subsetByColumn} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object 
#' @param colIndicator A \linkS4class{MultiAssayView} class object to be used for subsetting
#' @param drop logical indicates whether to coerce to lowest possible dimension after subsetting
subsetByColumn <- function(MultiAssayExperiment, colIndicator, drop = FALSE){
  if(is(colIndicator, "MultiAssayView") && getElement(colIndicator, "type") != "colnames"){
    stop("MultiAssayView class should be of colnames")
  } else {
    newMap <- getMap(colIndicator)
    subsetor <- lapply(colIndicator@keeps, function(x) unlist(x[,2]))
  }
  ##  mendoapply not working
  ## 	newSubset <- S4Vectors::mendoapply(subsetSample, Elist(MultiAssayExperiment), subsetor) 
  newSubset <- mapply(function(x, y, i){x[, y, drop = i]}, x = Elist(MultiAssayExperiment), y = subsetor, i = drop)
  newSubset <- Elist(newSubset)
  # Clone or replace method for slot??
  MultiAssayExperiment@sampleMap <- newMap
  MultiAssayExperiment@Elist <- newSubset
  MultiAssayExperiment@drops <- c(MultiAssayExperiment@drops, list(.convertList(colIndicator, "drops")))
  return(MultiAssayExperiment)
}
