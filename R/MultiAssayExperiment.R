### ==============================================
### MultiAssayExperiment constructor
### ----------------------------------------------

.generateMap <- function(mPheno, exlist) {
    samps <- lapply(exlist, colnames)
    listM <- lapply(seq_along(samps), function(i, x) {
    S4Vectors::DataFrame(assay = factor(names(x)[i]), colname = x[[i]])
  }, x = samps)
  full_map <- do.call(S4Vectors::rbind, listM)
  matches <- match(full_map$colname, rownames(mPheno))
  if (all(is.na(matches))) {
    stop("no way to map pData to ExperimentList")
  }
  primary <- rownames(mPheno)[matches]
  autoMap <- S4Vectors::cbind(full_map["assay"], DataFrame(primary),
                              full_map["colname"])
  if (any(is.na(autoMap$primary))) {
    notFound <- autoMap[is.na(autoMap$primary), ]
    warning("Data from rows:",
            sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
            "\ndropped due to missing phenotype data")
  }
  autoMap <- autoMap[!is.na(autoMap$primary), ]
  return(autoMap)
}

#' Create a MultiAssayExperiment object 
#'
#' This is the constructor function for the \link{MultiAssayExperiment-class}.
#' It combines multiple data elements from the different hierarchies of data
#' (study, experiments, and samples). It can create instances where neither
#' a \code{sampleMap} or a \code{pData} set is provided. Please see the
#' MultiAssayExperiment API documentation for more information by running the
#' \code{API} function.
#' 
#' @param ExperimentList A \code{list} or \link{ExperimentList} of all
#' combined experiments
#' @param pData A \code{\link[S4Vectors]{DataFrame}} or \code{data.frame} of
#' the phenotype data for all participants
#' @param sampleMap A \code{DataFrame} or \code{data.frame} of assay names,
#' sample identifiers, and colname samples
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
    function(experiments = ExperimentList(),
            pData = S4Vectors::DataFrame(),
            sampleMap = S4Vectors::DataFrame(),
            drops = list()) {
        if (inherits(experiments, "list"))
            experiments <- ExperimentList(experiments)
        else if (!inherits(experiments, "SimpleList"))
            stop("'experiments' must be a list or ExperimentList")
        if (!is(pData, "DataFrame"))
            pData <- S4Vectors::DataFrame(pData)
        if (!is(sampleMap, "DataFrame"))
            sampleMap <- S4Vectors::DataFrame(sampleMap)
        if (!all(c(length(sampleMap) == 0L,
                    length(pData) == 0L,
                    length(experiments) == 0L))) {
            if ((length(sampleMap) == 0L) && (length(pData) == 0L)) {
                allsamps <- unique(unlist(lapply(experiments, colnames)))
                pData <- S4Vectors::DataFrame(row.names = allsamps)
                sampleMap <- .generateMap(pData, experiments)
            } else if ((length(sampleMap) == 0L) && (length(pData) != 0L)) {
                sampleMap <- .generateMap(pData, experiments)
                validAssays <-
                    S4Vectors::split(
                        sampleMap[["colname"]], sampleMap[, "assay"])
                experiments <- Map(function(x, y) {
                    x[, y]
                }, experiments, validAssays)
                experiments <- ExperimentList(experiments)
            }
        }
        if (any(vapply(sampleMap, FUN = function(col) {
            !is.character(col)
        }, FUN.VALUE = logical(1L)))) {
            sampleMap[] <- lapply(sampleMap, as.character)
        }
        newMultiAssay <- new("MultiAssayExperiment",
                             ExperimentList = experiments,
                             pData = pData, 
                             sampleMap = sampleMap)
        return(newMultiAssay)
    }
