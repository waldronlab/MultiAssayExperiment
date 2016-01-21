#' Subset MultiAssayExperiment object
#' 
#' \code{subsetByColumn} returns a subsetted 
#' \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object 
#' @param colIndicator Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what rownames in the pData to select
#' for subsetting
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByColumn <- function(MultiAssayExperiment, colIndicator) {
  identifiers <- rownames(pData(MultiAssayExperiment))
  if(is.character(colIndicator)){
    logMatches <- identifiers %in% colIndicator
    if(!any(logMatches)){
      stop("No matching identifiers found")
    }
    selectors <- identifiers[identifiers %in% colIndicator]
  } else {
    selectors <- rownames(pData(MultiAssayExperiment))[colIndicator]
  }
  listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
  listMap <- listMap[order(names(MultiAssayExperiment))]
  listMap <- lapply(listMap, function(assay) {
    assay[which(as.vector(assay[, 1]) %in% selectors),]
  })
  newMap <- .convertList(listMap)
  columns <- lapply(listMap, function(x) {x[, 2, drop = TRUE]})
  newSubset <- mapply(function(x, j) {x[, j, drop = FALSE]},
                      x = Elist(MultiAssayExperiment), j = columns)
  newSubset <- Elist(newSubset)
  Elist(MultiAssayExperiment) <- newSubset
  sampleMap(MultiAssayExperiment) <- newMap
  return(MultiAssayExperiment)
} 
