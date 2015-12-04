.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.outersect <- function(x, y){
  c(setdiff(x, y), setdiff(y, x))
}

.separateMap <- function(object, ids){
	DFsampleMap <- S4Vectors::DataFrame(object@sampleMap)
	listDFsampleMap <- toListMap(DFsampleMap, "assayname")
	listDFsampleMap <- listDFsampleMap[order(names(object@elist))]
	loglistmatch <- lapply(listDFsampleMap, function(map) { map[,"master"] %in% ids })
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
#' @return A \code{\linkS4class{stage}} class object for subsequent subsetting
#' @export stage
stage <- function(MultiAssay, identifier, method = character(), ...){
  method <- match.arg(method, c("samples", "features", "assays"))
  if(method == "samples"){
    totalSamples <- samples(MultiAssay)
    if(!is.numeric(identifier) && !all(identifier %in% rownames(myMultiAssay@masterPheno))){
      iders <- intersect(identifier, rownames(MultiAssay@masterPheno))
      notUsed <- setdiff(identifier, rownames(MultiAssay@masterPheno))
      warning("Non-matched identifers will be dropped! : ", notUsed)
    } else {
      iders <- rownames(.subPheno(MultiAssay, identifier))
    }
    biMap <- .separateMap(MultiAssay, iders)
    newStage <- new("stage",
                    query = identifier,
                    keeps = biMap[["keeps"]],
                    drops = biMap[["drops"]],
                    type = "samples")
  } else if(method == "features"){
    totalFeatures <- features(MultiAssay)
    subsetor <- getHits(MultiAssay, identifier)
    newDrops <- .featMap(Map(function(x, y){.outersect(x, y)}, subsetor, totalFeatures))
    newKeeps <- .featMap(subsetor)
    newStage <- new("stage", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops,
                    type = "features")
  } else if(method == "assays"){
    if(!all(identifier %in% names(MultiAssay))){
      stop("Invalid experiment names!")
    }
    if(is.character(identifier)){
      newKeeps <- as.list(names(MultiAssay) %in% identifier)
    } else if(is.numeric(identifier)){
      newKeeps <- as.list(names(MultiAssay) %in% names(MultiAssay)[identifier])
    }
    names(newKeeps) <- names(MultiAssay)
    newDrops <- lapply(newKeeps, `!`) 
    newStage <- new("stage", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops, 
                    type = "assay")
  }
  return(newStage)
}
