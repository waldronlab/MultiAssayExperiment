.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.arrangeMap <- function(map, ordering){
	newOrd <- do.call(c, lapply(ordering, function(ord) { 
					 which(as.vector(map[, 1]) == ord) } ))
	return(newOrd)
}
.outersect <- function(x, y){
  c(setdiff(x, y), setdiff(y, x))
}

.separateMap <- function(object, ids){
	listDFsampleMap <- toListMap(object@sampleMap, "assayname")
	listDFsampleMap <- listDFsampleMap[order(names(object@Elist))]
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
      S4Vectors::DataFrame(feature = NA) # , assayname = names(x)[i])
    } else {
      S4Vectors::DataFrame(feature = x[[i]]) # , assayname = Rle(names(x)[i]))
    }
  }, x = featList)
  names(defeatmap) <- names(featList)
  return(defeatmap)
}

#' Stage by Samples, Features, or Assays
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}}
#' @param identifer Either a \code{character}, \code{numeric} or \code{logical} vector identifying targets 
#' @param method Prepare/Stage for subsetting using samples, features or assays.
#' @return A \code{\linkS4class{Stage}} class object for subsequent subsetting
#' @export Stage
Stage <- function(MultiAssay, identifier, method = character(), ...){
  method <- match.arg(method, c("samples", "features", "assays"))
  if(method == "samples"){
    totalSamples <- samples(MultiAssay)
    if(!is.numeric(identifier) && !all(identifier %in% rownames(myMultiAssay@masterPheno))){
      iders <- intersect(identifier, rownames(MultiAssay@masterPheno))
      notUsed <- setdiff(identifier, rownames(MultiAssay@masterPheno))
      warning("Nonmatched identifers were dropped! : ", notUsed)
    } else {
      iders <- rownames(.subPheno(MultiAssay, identifier))
    }
    biMap <- .separateMap(MultiAssay, iders)
    newStage <- new("Stage",
                    query = iders,
                    keeps = biMap[["keeps"]],
                    drops = biMap[["drops"]],
                    type = "samples")
  } else if(method == "features"){
    totalFeatures <- features(MultiAssay)
    subsetor <- getHits(MultiAssay, identifier)
    newDrops <- .featMap(Map(function(x, y){.outersect(x, y)}, subsetor, totalFeatures))
    newKeeps <- .featMap(subsetor)
    newStage <- new("Stage", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops,
                    type = "features")
  } else if(method == "assays"){
    if(is.logical(identifier)){
      if(length(identifier) == length(MultiAssay)){
      newKeeps <- identifier
    } else {
      stop("Provide a valid logical assay identifier of equal length!")
    }
    } else if(is.character(identifier)){
      if(all(identifier %in% names(MultiAssay))){
      newKeeps <- as.list(names(MultiAssay) %in% identifier)
    } else {
      stop("Provide a vector of valid experiment names!")
    }
    } else if(is.numeric(identifier)){
      if(all(identifier %in% 1:length(MultiAssay))){
      newKeeps <- as.list(names(MultiAssay) %in% names(MultiAssay)[identifier])
      } else {
        stop("Identifier out of bounds!")
      }
    }
    names(newKeeps) <- names(MultiAssay)
    newDrops <- lapply(newKeeps, `!`) 
    newStage <- new("Stage", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops, 
                    type = "assays")
  }
  return(newStage)
}
