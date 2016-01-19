#' Subset MultiAssayExperiment object
#' 
#' \code{subsetByColumn} returns a subsetted 
#' \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param MultiAssayExperiment A \code{\link{MultiAssayExperiment}} object 
#' @param colIndicator A \linkS4class{MultiAssayView} class object to be used 
#' for subsetting
#' @return A \code{\link{MultiAssayExperiment}} object 
subsetByColumn <- function(MultiAssayExperiment, colIndicator) {
    if (is.character(colIndicator)) {
        loglistmatch <- lapply(colnames(MultiAssayExperiment), 
                          function(charElem) {
                            charElem %in% colIndicator
                            })
        if (!any(unlist(loglistmatch))) {
            stop("No matches found")
        }
    } else {
        stop("colIndicator must be character")
    }
    newSubset <- mapply(function(x, i, j, drop) {
        x[, j, drop = FALSE]
    }, x = Elist(MultiAssayExperiment), j = loglistmatch)
    newSubset <- Elist(newSubset)
    listMap <- toListMap(sampleMap(MultiAssayExperiment), "assayname")
    logmapInd <- lapply(listMap, function(x) {
        x[, 2] %in% colIndicator
    })
    newMap <- mapply(function(x, y) {
        x[y, ]
    }, listMap, logmapInd)
    newMap <- Filter(function(x) {
        !isEmpty(x)
    }, newMap)
    newMap <- .convertList(newMap)
    sampleMap(MultiAssayExperiment) <- newMap
    Elist(MultiAssayExperiment) <- newSubset
    return(MultiAssayExperiment)
} 
