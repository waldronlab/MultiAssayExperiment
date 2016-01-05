.arrangeMap <- function(map, ordering){
	newOrd <- do.call(c, lapply(ordering, function(ord) { 
					 which(as.vector(map[, 1]) == ord) } ))
	return(newOrd)
}

.outersect <- function(x, y){
  c(setdiff(x, y), setdiff(y, x))
}

.separateMap <- function(object, ids){
	listDFsampleMap <- toListMap(sampleMap(object), "assayname")
	listDFsampleMap <- listDFsampleMap[order(names(Elist(object)))]
	loglistmatch <- lapply(listDFsampleMap, function(map) { as.vector(map[,"master"]) %in% ids })
	keeps <- Map(function(x, y) { x[y,] }, listDFsampleMap, loglistmatch)
	orderIndex <- lapply(keeps, .arrangeMap, ids)
	orderedKeeps <- Map(function(x, y) { x[y, ] }, keeps, orderIndex)
	duoMap <- list(keeps = orderedKeeps,
				   drops = Map(function(x, y) { x[y,] }, listDFsampleMap, lapply(loglistmatch, `!`)))
	return(duoMap)
}

.featMap <- function(featList){
  defeatmap <- lapply(seq_along(featList), FUN = function(i, x) {
    if(length(x[[i]])==0L){ 
      S4Vectors::DataFrame(feature = character(0)) # , assayname = names(x)[i])
    } else {
      S4Vectors::DataFrame(feature = x[[i]]) # , assayname = Rle(names(x)[i]))
    }
  }, x = featList)
  names(defeatmap) <- names(featList)
  return(defeatmap)
}

#' MultiAssayView operation for colnames, rownames, or assay
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}}
#' @param identifier Either a \code{character}, \code{numeric} or \code{logical} vector identifying targets 
#' @param method Prepare/View for subsetting using colnames, rownames or assays.
#' @param ... Additional arguments passed to the findOverlaps function
#' @return A \code{\linkS4class{MultiAssayView}} class object for subsequent subsetting
#' @export MultiAssayView
MultiAssayView <- function(MultiAssay, identifier, method = character(), ...){
  method <- match.arg(method, c("colnames", "rownames", "assays"))
  if(method == "colnames"){
    totalColnames <- colnames(MultiAssay)
    if(!is.numeric(identifier) && !all(identifier %in% rownames(masterPheno(myMultiAssay)))){
      iders <- intersect(identifier, rownames(masterPheno(MultiAssay)))
      notUsed <- setdiff(identifier, rownames(masterPheno(MultiAssay)))
      warning("Nonmatched identifers were dropped! : ", notUsed)
    } else {
      iders <- rownames(masterPheno(MultiAssay)[identifier, ])
    }
    biMap <- .separateMap(MultiAssay, iders)
    newMultiAssayView <- new("MultiAssayView",
                    query = iders,
                    keeps = biMap[["keeps"]],
                    drops = biMap[["drops"]],
                    type = "colnames")
  } else if(method == "rownames"){
    totalRownames <- rownames(MultiAssay)
    subsetor <- getHits(MultiAssay, identifier, ...)
    newDrops <- .featMap(Map(function(x, y){.outersect(x, y)}, subsetor, totalRownames))
    newKeeps <- .featMap(subsetor)
    newMultiAssayView <- new("MultiAssayView", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops,
                    type = "rownames")
  } else if(method == "assays"){
    if(is.logical(identifier)){
      if(length(identifier) == length(MultiAssay)){
        newKeeps <- as.list(identifier)
      } else {
        stop("Provide a valid logical assay identifier of equal length")
      }
    } else if(is.character(identifier)){
      if(all(identifier %in% names(MultiAssay))){
        newKeeps <- as.list(names(MultiAssay) %in% identifier)
      } else {
        stop("Provide a vector of valid experiment names")
      }
    } else if(is.numeric(identifier)){
      if(all(identifier %in% 1:length(MultiAssay))){
        newKeeps <- as.list(names(MultiAssay) %in% names(MultiAssay)[identifier])
      } else {
        stop("Identifier out of bounds")
      }
    }
    names(newKeeps) <- names(MultiAssay)
    newDrops <- lapply(newKeeps, `!`) 
    newMultiAssayView <- new("MultiAssayView", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops, 
                    type = "assays")
  }
  return(newMultiAssayView)
}
