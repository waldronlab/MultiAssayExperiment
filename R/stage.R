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
    if(is.character(identifier)){
      subsetor <- Map(function(x, y){x %in% y}, totalFeatures, identifier)
      newKeeps <- .featMap(Map(function(x, y) {x[y]}, totalFeatures, subsetor))
      newDrops <- .featMap(Map(function(x, y) {x[y]}, totalFeatures, lapply(subsetor, "!")))
      newStage <- new("stage", 
                      query = identifier, 
                      keeps = newKeeps,
                      drops = newDrops,
                      type = "features")
      return(newStage)
    } else if(is(identifer, "GRanges")){
      elist_classes <- sapply(MultiAssay@elist, class)
      logic_flag <- elist_classes %in% c("GRanges", "GRangesList", "RangedSummarizedExperiment")
      rangeBased <- Map(subset, MultiAssay@elist, logic_flag)
      rangeIdentifiers <- lapply(rangeBased, identify)      
      ## TODO: GRanges ---		findOverlaps(MultiAssay@elist, identifier, ...)@subjectHits
      ## TODO: GRangesList ---		lapply(rangedFeats, function(x) { findOverlaps(x, identifier, ...) } )
      ## TODO: RangedSummarizedExperiment --- findOverlaps(rowRanges(MultiAssay@elist), identifier, ...)@subjectHits
    }
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
