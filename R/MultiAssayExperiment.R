### ==============================================
### MultiAssayExperiment constructor
### ----------------------------------------------

.generateMap <- function(pData, experiments) {
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    matches <- match(colname, rownames(pData))
    if (length(matches) && all(is.na(matches)))
        stop("no way to map pData to ExperimentList")
    primary <- rownames(pData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)

    if (nrow(autoMap) && any(is.na(autoMap$primary))) {
        notFound <- autoMap[is.na(autoMap$primary), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap$primary), ]
    }
    autoMap
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
#' @param experiments A \code{list} or \link{ExperimentList} of all
#' combined experiments
#' @param pData A \code{\link[S4Vectors]{DataFrame}} or \code{data.frame} of
#' the phenotype data for all participants
#' @param sampleMap A \code{DataFrame} or \code{data.frame} of assay names,
#' sample identifiers, and colname samples
#' @param metadata An optional argument of "ANY" class (usually list) for
#' content describing the overall experiments.
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
            metadata = NULL,
            drops = list()) {
        if (inherits(experiments, "list"))
            experiments <- ExperimentList(experiments)
        else if (!inherits(experiments, "SimpleList"))
            stop("'experiments' must be a list or ExperimentList")
        if (!is(pData, "DataFrame"))
            pData <- S4Vectors::DataFrame(pData)
        if (!is(sampleMap, "DataFrame"))
            sampleMap <- S4Vectors::DataFrame(sampleMap)
        if (!all(c(ncol(sampleMap) == 0L,
                    ncol(pData) == 0L,
                    length(experiments) == 0L))) {
            if ((ncol(sampleMap) == 0L) && (ncol(pData) == 0L)) {
                allsamps <- unique(unlist(unname(colnames(experiments))))
                pData <- S4Vectors::DataFrame(row.names = allsamps)
                sampleMap <- .generateMap(pData, experiments)
            } else if ((ncol(sampleMap) == 0L) && (ncol(pData) != 0L)) {
                sampleMap <- .generateMap(pData, experiments)
                validAssays <-
                    S4Vectors::split(
                        sampleMap[["colname"]], sampleMap[, "assay"])
                experiments <- mapply(function(x, y) {
                    x[, y, drop = FALSE]
                }, experiments, validAssays, SIMPLIFY = FALSE)
                experiments <- ExperimentList(experiments)
            }
        }

        newMultiAssay <- new("MultiAssayExperiment",
                             ExperimentList = experiments,
                             pData = pData,
                             sampleMap = sampleMap,
                             metadata = metadata)
        return(newMultiAssay)
    }
