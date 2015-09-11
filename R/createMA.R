#' Create a MultiAssayExperiment object 
#' \code{createMA} returns a \code{\linkS4class{MultiAssayExperiment}} object 
#'
#' This function combines multiple data sources specific to one disease by matching samples. 
#' 
#' @param masterpheno A data frame of the phenotype data for all participants.
#' @param objlist A list of all combined experiments
#' @param drop Logical (default FALSE) parameter for dropping samples with unmatched phenotype data.   
#' @return A \code{MultiAssayExperiment} data object that stores experiment and phenotype data.
createMA <- function(masterpheno, objlist, drop=FALSE, samplemaps=NULL){
    ## samplemaps will be maps that rename samples in object list to names used in masterpheno.
	if(!is(masterpheno, "data.frame")){
		stop("masterpheno should be a data frame of metadata for all samples")
	}
	if(!is(objlist, "list")){
		stop("objlist should be a named list of data objects")
	}

	if(!is.null(samplemaps)){
		if(!any(c(is(samplemaps, "data.frame"), is(samplemaps, "list")))){
			stop("samplemaps should either be a list or a data frame!")
		} else if(is(samplemaps, "list")){
			# if(!all(grepl("[A-Z]{4}.[0-9]{2}.[0-9]{4}", unlist(samplemaps))))
			if(all(unique(as.character(unname(unlist(sapply(samplemaps, "[", 1)))))) %in% rownames(masterpheno)){
				message("All unique IDs matched to the masterpheno data frame!")
			} else {
				warning("All IDs not present in the masterpheno!")
			}
		} else if(is(samplemaps, "data.frame")){
			if(all(unique(as.character(samplemaps[,1])) %in% rownames(masterpheno))){
				message("All unique IDs matched to the masterpheno data frame!")
			} else {
				message("All IDs not present in the masterpheno!")
			}
			#' 			availableIDs <-	samplemaps[which(!is.na(samplemaps))]
			#' 			if(!all(grepl("[A-Z]{4}.[0-9]{2}.[0-9]{4}", availableIDs))){
			#' 				stop("All available sample names in the data frame must have a specific alphanumeric format!")}
		}
	} else { warning("No sample maps provided!") } 

if(is(samplemaps, "list")){
inpheno <- lapply(samplemaps, function(exptmap) rownames(exptmap) %in% rownames(masterpheno))
}else if(is(samplemaps, "data.frame")){
inpheno <- apply(samplemaps, 2, FUN = function(expt) { expt %in% rownames(masterpheno) } )
}	

    ##Sample names checking:
    has.pheno <- lapply(objlist, function(x) colnames(x) %in% rownames(masterpheno))
    if(!drop){
	errmsg <- paste("Missing the following number of masterpheno entries for each data type: ",
			paste(names(objlist), ":", sapply(has.pheno, function(x) sum(!x)), collapse=", "),
			". Set drop=TRUE to drop these observations, or add samples to masterpheno.")
	stop(errormsg)
    }else{
	message("Dropping the following samples:")
	for (i in 1:length(objlist)){
	    if(all(has.pheno[[i]])) next
	    message(paste(names(objlist)[i], ":", collapse=""))
	    message(paste(colnames(objlist[[i]])[!has.pheno[[i]]], collapse=" "))
	    message("\n ")
	    objlist[[i]] <- objlist[[i]][, has.pheno[[i]]]
	}

    }
    exptlist <- lapply(1:length(objlist), function(i) 
	new("expt", serType="in-memory", assayPath="", tag=names(objlist)[i]))
    hub <- new("eHub", hub=exptlist, masterSampleData=masterpheno)
    res <- new("MultiAssayExperiment", basehub=hub, elist=objlist)
}
