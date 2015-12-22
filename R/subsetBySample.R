#' Subset MultiAssayExperiment object
#' \code{subsetBySample} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssay A \code{\link{MultiAssayExperiment}} object 
#' @param identifier A \linkS4class{Stage} class object to be used for subsetting
#' @export subsetBySample
subsetBySample <- function(MultiAssay, identifier){
  if(is(identifier, "Stage") && getElement(identifier, "type") != "samples"){
    stop("Stage class should be of samples!")
  } else {
    newMap <- getMap(identifier)
    subsetor <- lapply(identifier@keeps, function(x) unlist(x[,2]))
  }
  ##  mendoapply not working
  ## 	newSubset <- S4Vectors::mendoapply(subsetSample, MultiAssay@Elist, subsetor) 
  newSubset <- mapply(function(x, y){ x[ ,y , drop = FALSE]}, MultiAssay@Elist, subsetor)
  newSubset <- Elist(newSubset)
  # Clone or replace method for slot??
  MultiAssay@sampleMap <- newMap
  MultiAssay@Elist <- newSubset
  MultiAssay@drops <- c(MultiAssay@drops, list(.convertList(identifier, "drops")))
  return(MultiAssay)
}
