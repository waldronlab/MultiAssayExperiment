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
#' @param explist A \code{list} of all combined experiments
#' @param masterPheno A \code{data.frame} of the phenotype data for all participants.
#' @param sampleMap A \code{data.frame} of sample identifiers, assay samples, and assay names.
#' @param drop Logical (default FALSE) parameter for dropping samples with unmatched phenotype data.   
#' @return A \code{MultiAssayExperiment} data object that stores experiment and phenotype data.
#' @export MultiAssayExperiment
MultiAssayExperiment <- function(explist = list(), masterPheno = data.frame(), sampleMap = data.frame(), drops = list()){
	if(!all(c(length(sampleMap) == 0L, length(masterPheno) == 0L, length(explist) == 0L))){
		if((length(sampleMap) == 0L) & (length(masterPheno) == 0L)){
			allsamps <- unique(unlist(lapply(explist, samples)))
			masterPheno <- data.frame(pheno1 = rep(NA, length(allsamps)), row.names = allsamps, stringsAsFactors = FALSE)
			sampleMap <- .generateMap(masterPheno, explist)
		} else if((length(sampleMap) == 0L) & !(length(masterPheno) == 0L)){
			warning("sampleMap not provided! Map will be created from data provided.")
			sampleMap <- .generateMap(masterPheno, explist)
			validAssays <- split(sampleMap[["assay"]], sampleMap$assayname)
			explist <- Map(subsetSample, explist, validAssays) 
		}
		if(is(explist, "list")){
			explist <- S4Vectors::SimpleList(explist)
		}
	} else {
		explist <- S4Vectors::SimpleList(explist)
	}
	newMultiAssay <- new("MultiAssayExperiment",
						 elist = explist,
						 masterPheno = masterPheno,
						 sampleMap = sampleMap)
	return(newMultiAssay)
}

#if(!drop){
#	errmsg <- paste("Missing the following number of masterPheno entries for each data type: ",
#					paste(names(objlist), ":", sapply(has.pheno, function(x) sum(!x)), collapse=", "),
#					". Set drop=TRUE to drop these observations, or add samples to masterPheno.")
#	stop(errormsg)
#}else{
#	message("Dropping the following samples:")
#	for (i in 1:length(objlist)){
#		if(all(has.pheno[[i]])) next
#		message(paste(names(objlist)[i], ":", collapse=""))
#		message(paste(colnames(objlist[[i]])[!has.pheno[[i]]], collapse=" "))
#		message("\n ")
#		objlist[[i]] <- objlist[[i]][, has.pheno[[i]]]
#	}
#}
#	exptlist <- lapply(1:length(objlist), function(i) 
#					   new("expt", serType="in-memory", assayPath="", tag=names(objlist)[i]))
#	hub <- new("eHub", hub=exptlist, masterSampleData=masterPheno)
#	res <- new("MultiAssayExperiment", basehub=hub, elist=objlist)
#}
