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

.createNames <- function(object){
  if(class(object) %in% c("RaggedRangedAssay", "GRangesList")){
    for(i in seq_along(object)){
      names(object[[i]]) <- 1:length(object[[i]])
    }
  } else if(obj_cl == "RangedSummarizedExperiment"){
    names(object) <- 1:length(object)
  }
  return(object)
}

.PrepElements <- function(object){
  if(is.null(rownames(object))){
    object <- .createNames(object)
  }
  if(class(object) == "GRangesList"){
    object <- RaggedRangedAssay(object)
  } else { object } 
  return(object)
}

.generateMap <- function(mPheno, exlist){
  samps <- lapply(exlist, colnames)
  # listM <- lapply(seq_along(samps), function(i, x) {data.frame(assay = x[[i]], assayname = names(x)[i], row.names = NULL, stringsAsFactors = FALSE)}, x = samps)
  listM <- lapply(seq_along(samps), function(i, x) {S4Vectors::DataFrame(assay = x[[i]], assayname = names(x)[i])}, x = samps)
  full_map <- do.call(S4Vectors::rbind, listM)
  # full_map <- do.call(rbind, listM)
  master <- rownames(mPheno)[match(full_map$assay, rownames(mPheno))]
  # autoMap <- cbind(master, full_map)
  autoMap <- S4Vectors::cbind(DataFrame(master), full_map)
  if(any(is.na(autoMap$master))){
    notFound <- autoMap[is.na(autoMap$master),]
    warning("Data from rows:", sprintf("\n %s - %s", notFound[, 2], notFound[, 3]), "\ndropped due to missing phenotype data!")
  }
  autoMap <- autoMap[!is.na(autoMap$master),]
  return(autoMap)
}


#' Create a MultiAssayExperiment object 
#'
#' This function combines multiple data sources specific to one disease by matching samples. 
#' 
#' @param Elist A \code{list} of all combined experiments
#' @param masterPheno A \code{\link[S4Vectors]{DataFrame-class}} of the phenotype data for all participants.
#' @param sampleMap A \code{data.frame} of sample identifiers, assay samples, and assay names.
#' @param drops A \code{list} of unmatched information (included after subsetting)   
#' @return A \code{MultiAssayExperiment} data object that stores experiment and phenotype data.
#' @export MultiAssayExperiment
MultiAssayExperiment <- function(Elist = list(), masterPheno = S4Vectors::DataFrame(), sampleMap = S4Vectors::DataFrame(), drops = list()){
  Elist <- lapply(Elist, function(x) {.PrepElements(x)})
  if(!all(c(length(sampleMap) == 0L, length(masterPheno) == 0L, length(Elist) == 0L))){
    if((length(sampleMap) == 0L) & (length(masterPheno) == 0L)){
      allsamps <- unique(unlist(lapply(Elist, colnames)))
      masterPheno <- S4Vectors::DataFrame(pheno1 = rep(NA, length(allsamps)), row.names = allsamps)
      sampleMap <- .generateMap(masterPheno, Elist)
    } else if((length(sampleMap) == 0L) & !(length(masterPheno) == 0L)){
      warning("sampleMap not provided! Map will be created from data provided...")
      sampleMap <- .generateMap(masterPheno, Elist)
      validAssays <- split(sampleMap[["assay"]], sampleMap$assayname)
      Elist <- Map(function(x, y) { x[, y]}, Elist, validAssays) 
    }
  }
  if(!is(masterPheno, "DataFrame")){masterPheno <- S4Vectors::DataFrame(masterPheno)}
  if(!is(sampleMap, "DataFrame")){sampleMap <- S4Vectors::DataFrame(sampleMap)}
  newMultiAssay <- new("MultiAssayExperiment",
                       Elist = Elist(Elist),
                       masterPheno = masterPheno,
                       sampleMap = sampleMap)
  return(newMultiAssay)
}
