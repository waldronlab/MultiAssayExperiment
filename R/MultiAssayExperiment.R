.fillMap <- function(mapsection, phenonames){
	if(all(phenonames %in% unique(mapsection[,1]))){
		return(mapsection)
	} else {
		newrows <- phenonames[!(phenonames %in% unique(mapsection[,1]))]
		newmapsec <- rbind(mapsection,
						   matrix(c(newrows, rep(NA, length(newrows))),
								  ncol = 2,
								  dimnames = list(NULL, colnames(mapsection))))
		return(newmapsec)
	}
}

.ldmap <- function(listmap){
	dfmap <- lapply(seq_along(listmap), FUN = function(i, x) {
					data.frame(master = x[[i]][, 1],
							   assay = x[[i]][, 2],
							   assayname = names(x)[i],
							   row.names = NULL,
							   stringsAsFactors = FALSE) }, x = listmap)
	return(do.call(rbind, dfmap))
}

.FixElemNames <- function(object){
  obj_cl <- class(object)
  if(obj_cl == "GRangesList"){
    if(is.null(names(object[[1]]))){
      object <- createNames(object)
    } 
  } else { object } 
  return(object)
}

.generateMap <- function(mPheno, exlist){
	samps <- lapply(exlist, samples)
	listM <- lapply(seq_along(samps), function(i, x) {data.frame(assay = x[[i]], assayname = names(x)[i], row.names = NULL, stringsAsFactors = FALSE) }, x = samps)
	full_map <- do.call(rbind, listM)
	master <- rownames(mPheno)[match(full_map$assay, rownames(mPheno))]
	autoMap <- cbind(master, full_map)
	if(any(is.na(autoMap$master))){
		notFound <- autoMap[is.na(autoMap$master),]
		warning("Data from rows:", sprintf("\n %s - %s", notFound[, 2], notFound[, 3]), "\ndropped due to missing phenotype data!")
	}
	autoMap <- autoMap[!is.na(autoMap$master),]
	return(autoMap)
}

#' Create a MultiAssayExperiment object 
#' \code{MultiAssayExperiment} returns a \code{\linkS4class{MultiAssayExperiment}} object 
#'
#' This function combines multiple data sources specific to one disease by matching samples. 
#' 
#' @param elist A \code{list} of all combined experiments
#' @param masterPheno A \code{data.frame} of the phenotype data for all participants.
#' @param sampleMap A \code{data.frame} of sample identifiers, assay samples, and assay names.
#' @param drops A \code{list} of unmatched information (included after subsetting)   
#' @return A \code{MultiAssayExperiment} data object that stores experiment and phenotype data.
#' @export MultiAssayExperiment
MultiAssayExperiment <- function(elist = list(), masterPheno = data.frame(), sampleMap = data.frame(), drops = list()){
  elist <- lapply(elist, .FixElemNames)
	if(!all(c(length(sampleMap) == 0L, length(masterPheno) == 0L, length(elist) == 0L))){
		if((length(sampleMap) == 0L) & (length(masterPheno) == 0L)){
			allsamps <- unique(unlist(lapply(elist, samples)))
			masterPheno <- data.frame(pheno1 = rep(NA, length(allsamps)), row.names = allsamps, stringsAsFactors = FALSE)
			sampleMap <- .generateMap(masterPheno, elist)
		} else if((length(sampleMap) == 0L) & !(length(masterPheno) == 0L)){
			warning("sampleMap not provided! Map will be created from data provided.")
			sampleMap <- .generateMap(masterPheno, elist)
			validAssays <- split(sampleMap[["assay"]], sampleMap$assayname)
			elist <- Map(subsetSample, elist, validAssays) 
		}
	} else {
		newMultiAssay <- new("MultiAssayExperiment", elist = elist(elist), masterPheno = masterPheno, sampleMap = sampleMap)
	}
	newMultiAssay <- new("MultiAssayExperiment",
						 elist = elist(elist),
						 masterPheno = masterPheno,
						 sampleMap = sampleMap)
	return(newMultiAssay)
}

