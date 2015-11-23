.subPheno <- function(object, j){
	return(object@masterPheno[j, ])
}

.separateMap <- function(object, ids){
	DFsampleMap <- S4Vectors::DataFrame(object@sampleMap)
	listDFsampleMap <- toListMap(DFsampleMap, "assayname")
	browser()
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
					data.frame(feature = x[[i]],
							   assayname = names(x)[i],
							   row.names = NULL,
							   stringsAsFactors = FALSE) }, x = featList)
	return(do.call(rbind, defeatmap))
}

#' Stage by Samples, Features, or Assays
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}}
#' @param identifer Either a \code{character}, \code{numeric} or \code{logical} vector identifying targets 
#' @param by Stage for subsetting by either samples, features or assays.
#' @return A \code{\linkS4class{stage}} class object for subsequent subsetting
#' @export stage
stage <- function(MultiAssay, identifier, by = NULL, ...){
  by <- tolower(gsub("s$", "", by, ignore.case = TRUE))
  if(by == "sample"){
    totalSamples <- samples(MultiAssay)
    if(is.numeric(identifier)){
      iders <- rownames(.subPheno(MultiAssay, identifier))
    } else {
      iders <- .subPheno(MultiAssay, identifier)
    }
    biMap <- .separateMap(MultiAssay, iders)
    newStage <- new("stage",
                    query = iders,
                    keeps = biMap[["keeps"]],
                    drops = biMap[["drops"]],
                    type = "samples")
    return(newStage)
  } else if(by == "feature"){
    totalFeatures <- features(MultiAssay)
    newKeeps <- stage(MultiAssay, identifier)
    subsetor <- Map(function(x, y){x %in% y}, totalFeatures, identifier)
    newDrops <- .featMap(Map(function(x, y) {x[y]}, totalFeatures, lapply(subsetor, "!")))
    newStage <- new("stage", 
                    query = identifier, 
                    keeps = newKeeps,
                    drops = newDrops,
                    type = "features")
    return(newStage)
  } else if(by == "assay"){
    if(is.character(identifier)){
      newKeeps <- names(MultiAssay@elist)[identifier]
      newDrops <- !(names(MultiAssay@elist) %in% newKeeps)
      newStage <- new("stage", 
                      query = identifier, 
                      keeps = newKeeps,
                      drops = newDrops, 
                      type = "assay")
      return(newStage)
    } else {
      stop("Please specify assays by their name!")
    }
  }
}
