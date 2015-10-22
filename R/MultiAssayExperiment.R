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
MultiAssayExperiment <- function(explist = list(), masterPheno = data.frame(), sampleMap = list(), fill = FALSE, drop=FALSE){
	if(length(sampleMap) == 0L){
		##
		## TODO: create auto sampleMap function here
		## 
	}

	if(is(explist, "list")){
		explist <- S4Vectors::SimpleList(explist)
	}
if(fill){
##
## TODO: add columns where missing data present
##
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
