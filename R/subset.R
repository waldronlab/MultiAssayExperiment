.convertList <- function(object, slot=NULL, type = "colnames"){
  if(!is.null(slot)){
    listmap <- getElement(object, slot)
    type <- object@type
  } else {
    listmap <- object
  }
  if(is(listmap, "MultiAssayView")){stop("Provide a slot name for the MultiAssayView class")}
  DFmap <- lapply(seq_along(listmap), FUN = function(i, x){
    if(type == "colnames"){
      if(isEmpty(x[i])){
        S4Vectors::DataFrame(master = Rle(NA), assay = NA, assayname = Rle(names(x)[i]))
      } else {
        S4Vectors::DataFrame(master = Rle(x[[i]][, 1]),
                             assay = x[[i]][, 2],
                             assayname = Rle(names(x)[i]))
      }
    } else if(type == "rownames"){
      if(isEmpty(x[i])){
        S4Vectors::DataFrame(feature = NA, assayname = Rle(names(x)[i]))
      } else {
        S4Vectors::DataFrame(feature = x[[i]][, 1], assayname = Rle(names(x)[i]))
      }
    } else if(type == "assays"){
      S4Vectors::DataFrame(value = x[[i]], assayname = names(x)[i])
    }
  }, x = listmap)
  newMap <- do.call(S4Vectors::rbind, DFmap)
  newMap <- newMap[!is.na(newMap[, 1]),]
  return(newMap)
}

#' Subset MultiAssayExperiment object
#' \code{subset} returns a subsetted \code{\linkS4class{MultiAssayExperiment}} object
#' 
#' @param x A \code{\linkS4class{MultiAssayExperiment}} object
#' @param indicator A \code{logical} or \code{character} vector or \code{MultiAssayView} class object to use for subsetting
#' @param method A \code{character} vector of length one designating to subset either by colnames, rownames, or assays
#' @param drop logical (default FALSE) whether to coerce lowest possible dimension after subsetting 
#' @param ... Additional arguments to pass to SubsetByOverlaps when subsetting by rownames
#' @describeIn MultiAssayExperiment
#' export subset
setMethod("subset", "MultiAssayExperiment", function(x, indicator, method = NULL, drop = FALSE, ...){
  if(is(indicator, "MultiAssayView")){
    method <- getElement(indicator, "type") 
  } else if(is.null(method)){
    stop("Indicate a subset method")
  } else {
    method <- match.arg(method, c("colnames", "rownames", "assays"))
  }
  if(method == "colnames"){
    return(subsetByColumn(MultiAssayExperiment = x, colIndicator = indicator, drop))
  } else if(method == "rownames"){
    return(subsetByRow(MultiAssayExperiment = x, rowIndicator = indicator, drop, ...))
  } else if(method == "assays"){
    return(subsetByAssay(MultiAssayExperiment = x, assayIndicator = indicator, drop))
  } 
})
