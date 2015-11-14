#' Stage by Samples, Features, or Assays
#' 
#' @param MultiAssay A \code{\linkS4class{MultiAssayExperiment}}
#' @param identifer Either a \code{character}, \code{numeric} or \code{logical} vector identifying targets 
#' @return A \code{\linkS4class{stage}} class object for subsequent subsetting
stage <- function(MultiAssay, identifier, by = NULL){
	if(by == "samples"){
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
	} else if(by == "features"){
		totalFeatures <- features(MultiAssay)
		charFeats <- totalFeatures[which(sapply(totalFeatures, class) == "character")]
		rangeFeats <- totalFeatures[!(names(MultiAssay) %in% names(charFeats))]
		if(is.character(identifier)){
			newKeeps <- lapply(charFeats, function(feats) { intersect(identifier, feats) })
			charDrops <- lapply(charFeats, function(feast) { sort(c(setdiff(identifier, feats), setdiff(feats, identifier))) })
			newDrops <- c(charDrops, rangeFeats)
			newStage <- new("stage", 
							query = identifier, 
							keeps = newKeeps
							drops = newDrops
							type = "features")
			return(newStage)
		} else if(is(identifer, "GRanges")){
			## TODO:			overlapsAny(rangeFeats, identifer)
		} else if(is(identifier, "GRangesList")) { 
			## TODO:		lapply(identifer, function(x) { overlapsAny(rangeFeats, x) } )
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
