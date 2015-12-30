#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param identifier A \linkS4class{MultiAssayView} class object to be used for subsetting
subsetBySample <- function(MultiAssay, identifier, drop = FALSE){
  if(is(identifier, "MultiAssayView") && getElement(identifier, "type") != "colnames"){
    stop("MultiAssayView class should be of colnames")
  } else {
    newMap <- getMap(identifier)
    subsetor <- lapply(identifier@keeps, function(x) unlist(x[,2]))
  }
  ##  mendoapply not working
  ## 	newSubset <- S4Vectors::mendoapply(subsetSample, Elist(MultiAssay), subsetor) 
  newSubset <- mapply(function(x, y, i){x[, y, drop = i]}, x = Elist(MultiAssay), y = subsetor, i = drop)
  newSubset <- Elist(newSubset)
  # Clone or replace method for slot??
  MultiAssay@sampleMap <- newMap
  MultiAssay@Elist <- newSubset
  MultiAssay@drops <- c(MultiAssay@drops, list(.convertList(identifier, "drops")))
  return(MultiAssay)
}
