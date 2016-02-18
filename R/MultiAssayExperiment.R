### ==============================================
### MultiAssayExperiment constructor
### ----------------------------------------------

.generateMap <- function(mPheno, exlist) {
  samps <- lapply(exlist, colnames)
  listM <- lapply(seq_along(samps), function(i, x) {
    S4Vectors::DataFrame(assay = x[[i]], assayname = Rle(names(x)[i]))
  }, x = samps)
  full_map <- do.call(S4Vectors::rbind, listM)
  master <- Rle(rownames(mPheno)[match(full_map$assay, rownames(mPheno))])
  autoMap <- S4Vectors::cbind(DataFrame(master), full_map)
  if (any(is.na(autoMap$master))) {
    notFound <- autoMap[is.na(autoMap$master), ]
    warning("Data from rows:",
            sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
            "\ndropped due to missing phenotype data")
  }
  autoMap <- autoMap[!is.na(autoMap$master), ]
  return(autoMap)
}

#' Create a MultiAssayExperiment object 
#'
#' This function combines multiple data elements from the different hierarchies
#' of data (study, experiments, and samples)
#' 
#' @param Elist A \code{list} of all combined experiments
#' @param pData A \code{\link[S4Vectors]{DataFrame-class}} of the phenotype
#' data for all participants
#' @param sampleMap A \code{DataFrame} of sample identifiers, assay samples,
#' and assay names
#' @param drops A \code{list} of unmatched information
#' (included after subsetting)   
#' @return A \code{MultiAssayExperiment} data object that stores experiment
#' and phenotype data
#' 
#' @example inst/scripts/MultiAssayExperiment-Ex.R
#' 
#' @export MultiAssayExperiment
#' @seealso MultiAssayExperiment-class
MultiAssayExperiment <-
  function(Elist = list(),
           pData = S4Vectors::DataFrame(),
           sampleMap = S4Vectors::DataFrame(),
           drops = list()) {
    newElist <- Elist(Elist)
    if (!all(c(length(sampleMap) == 0L,
               length(pData) == 0L,
               length(Elist) == 0L))) {
      if ((length(sampleMap) == 0L) && (length(pData) == 0L)) {
        allsamps <- unique(unlist(lapply(Elist, colnames)))
        pData <- S4Vectors::DataFrame(
          pheno1 = rep(NA, length(allsamps)),
          row.names = allsamps)
        sampleMap <- .generateMap(pData, newElist)
      } else if ((length(sampleMap) == 0L) && !(length(pData) == 0L)) {
        warning("sampleMap not provided, map will be generated")
        sampleMap <- .generateMap(pData, Elist)
        validAssays <-
          S4Vectors::split(sampleMap[["assay"]], sampleMap[, "assayname"])
        newElist <- Map(function(x, y) {
          x[, y]
        }, newElist, validAssays)
        newElist <- Elist(newElist)
      }
    }
    if (!is(pData, "DataFrame")) {
      pData <- S4Vectors::DataFrame(pData)
    }
    if (!is(sampleMap, "DataFrame")) {
      sampleMap <- S4Vectors::DataFrame(sampleMap)
    }
    newMultiAssay <- new("MultiAssayExperiment",
                         Elist = newElist,
                         pData = pData, 
                         sampleMap = sampleMap)
    return(newMultiAssay)
  }
