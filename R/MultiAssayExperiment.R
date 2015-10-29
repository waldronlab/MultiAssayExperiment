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

.generateMap <- function(masterPheno, exlist){
	samps <- lapply(exlist, samples)
	masterlist <- lapply(samps, FUN = function(x) { union(x, rownames(masterPheno)) })
	genMap <- Map(function(x, y) { data.frame(matrix(cbind(x, ifelse(x %in% y, x, NA)),
														nrow = length(x),
														dimnames = list(seq_along(x), c("master", "assay"))), stringsAsFactors = FALSE)}, 
				  masterlist, samps)
}

#' Create a MultiAssayExperiment object 
#' \code{MultiAssayExperiment} returns a \code{\linkS4class{MultiAssayExperiment}} object 
#'
#' This function combines multiple data sources specific to one disease by matching samples. 
#' 
#' @param explist A list of all combined experiments
#' @param masterPheno A data frame of the phenotype data for all participants.
#' @param sampleMap A list object of identifiers and sample names for each experiment
#' @param drop Logical (default FALSE) parameter for dropping samples with unmatched phenotype data.   
#' @return A \code{MultiAssayExperiment} data object that stores experiment and phenotype data.
#' @export MultiAssayExperiment
MultiAssayExperiment <- function(explist = list(), masterPheno = data.frame(), sampleMap = list(), drop=FALSE){
	if((length(sampleMap) == 0L) & (length(masterPheno) == 0L)){
		allsamps <- lapply(explist, samples)
		sampleMap <- lapply(allsamps, function(map) { data.frame(master = map, stringsAsFactors = FALSE) })
	} else if((length(sampleMap) == 0L) & !(length(masterPheno) == 0L)){
		warning("sampleMap not provided! Map will be created from data provided.")
		sampleMap <- .generateMap(masterPheno, explist)
		notFound <- lapply(sampleMap, FUN = function(x) { x[!(x[,1] %in% rownames(masterPheno)), 1] } )
		notFound <- Reduce(union, notFound)
		masterPheno <- rbind(masterPheno, matrix(NA,
												 nrow = length(notFound), ncol = length(masterPheno),
												 dimnames = list(notFound, names(masterPheno))
												 ))
	}
	sampleMap <- lapply(sampleMap, .fillMap, rownames(masterPheno))
	if(is(explist, "list")){
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
