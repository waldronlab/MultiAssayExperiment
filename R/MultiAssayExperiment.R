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

.harmonize <- function(experiments, pData, sampleMap) {
    harmony <- character()
    ## experiment and sampleMap assays need to agree
    assay <- intersect(names(experiments), levels(sampleMap[["assay"]]))
    keep_sampleMap_assay <- sampleMap[["assay"]] %in% assay
    if (!all(keep_sampleMap_assay)) {
        sampleMap <- sampleMap[keep_sampleMap_assay, , drop=FALSE]
        sampleMap[["assay"]] <- factor(sampleMap[["assay"]], levels=assay)
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_assay),
                  "sampleMap rows not in names(experiments)"))
    }

    ## experiment colnames and sampleMap colname need to agree
    grp <- sampleMap$assay
    colnm = split(sampleMap$colname, grp)
    keep = Map(intersect, colnm, colnames(experiments)[names(colnm)])
    keep_sampleMap_colname = logical(nrow(sampleMap))
    split(keep_sampleMap_colname, grp) <- Map("%in%", colnm, keep)
    if (!all(keep_sampleMap_colname)) {
        sampleMap <- sampleMap[keep_sampleMap_colname, , drop=FALSE]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_colname),
                  "sampleMap rows with 'colname' not in colnames of experiments"))
    }

    ## primary and sampleMap primary need to agree
    primary <- intersect(rownames(pData), sampleMap[["primary"]])
    keep_sampleMap_primary <- sampleMap[["primary"]] %in% primary
    if (!all(keep_sampleMap_primary)) {
        sampleMap <- sampleMap[keep_sampleMap_primary,, drop=FALSE]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_primary),
                  "sampleMap rows with 'primary' not in pData"))
    }

    ## update objects
    assay <- intersect(names(experiments), levels(sampleMap[["assay"]]))
    experiments_columns <- split(sampleMap[["colname"]], sampleMap[["assay"]])
    primary <- intersect(rownames(pData), sampleMap[["primary"]])
    keep_pData <- rownames(pData) %in% primary
    if (!all(keep_pData)){
        pData=pData[keep_pData,]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_pData),
                  "pData rownames not in sampleMap 'primary'"))
    }

    experiments <- ExperimentList(Map(function(x, idx) {
        x[, colnames(x) %in% idx, drop=FALSE]
    }, experiments[assay], experiments_columns[assay]))

    if (length(harmony))
        message("harmonizing input:\n  ", paste(harmony, collapse="\n  "))
    list(experiments=experiments, sampleMap=sampleMap, pData=pData)
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

        if (missing(experiments))
            experiments = ExperimentList()
        else
            experiments <- ExperimentList(experiments)


        if (missing(pData)){
            allsamps <- unique(unlist(unname(colnames(experiments))))
            pData <- S4Vectors::DataFrame(row.names = allsamps)
        } else if (!is(pData, "DataFrame"))
            pData <- S4Vectors::DataFrame(pData)


        if (missing(sampleMap)){
            sampleMap <- .generateMap(pData, experiments)
        } else {
            sampleMap <- S4Vectors::DataFrame(sampleMap)
            if (!all(c("assay", "primary", "colname") %in% colnames(sampleMap)))
                stop("'sampleMap' does not have required columns")
            if (!is.factor(sampleMap[["assay"]]))
                sampleMap[["assay"]] <- factor(sampleMap[["assay"]])
            if (!is.character(sampleMap[["primary"]])) {
                warning("sampleMap[['primary']] coerced to character()")
                sampleMap[["primary"]] <- as.character(sampleMap[["primary"]])
            }
            if (!is.character(sampleMap[["colname"]])) {
                warning("sampleMap[['colname']] coerced to character()")
                sampleMap[["colname"]] <- as.character(sampleMap[["colname"]])
            }
        }

        bliss <- .harmonize(experiments, pData, sampleMap)

        ## validAssays <- S4Vectors::split(
        ##     sampleMap[["colname"]], sampleMap[, "assay"])
        ## experiments <- ExperimentList(Map(function(x, y) {
        ##     x[, y]
        ## }, experiments, validAssays))


        newMultiAssay <- new("MultiAssayExperiment",
                             ExperimentList = bliss[["experiments"]],
                             pData = bliss[["pData"]],
                             sampleMap = bliss[["sampleMap"]],
                             metadata = metadata)
        return(newMultiAssay)
    }
