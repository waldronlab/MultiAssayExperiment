#' @include MultiAssayExperiment-subset.R
NULL

.warnnomatch <- function(type, values) {
    warning("'", type, "': ",
        paste(S4Vectors:::selectSome(values), collapse = ", "),
            "\n  could not be matched")
}

.sampleMapFromData <- function(colData, experiments) {
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    matches <- match(colname, rownames(colData))
    if (length(matches) && all(is.na(matches)))
        stop("no way to map colData to ExperimentList")
    primary <- rownames(colData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)
    missingPrimary <- is.na(autoMap[["primary"]])
    if (nrow(autoMap) && any(missingPrimary)) {
        notFound <- autoMap[missingPrimary, ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!missingPrimary, ]
    }
    autoMap
}

.sampleMapFromExisting <- function(colData, sampleMap, colnames) {
    colmatches <- match(colnames, sampleMap[["colname"]])
    notFounds <- which(is.na(colmatches))
    if (!length(colmatches) || all(is.na(colmatches))) {
        pmatch <- match(rownames(colData), colnames)
        noPrim <- which(is.na(pmatch))
        pnames <- colnames[na.omit(pmatch)]
        if (!length(pnames) || all(is.na(pnames)))
            stop("No 'colnames' in experiments could be matched:\n  ",
                 paste(S4Vectors:::selectSome(colnames), collapse = ", "))
        else if (length(noPrim))
            .warnnomatch("primary", colnames[noPrim])
        DataFrame(primary = pnames, colname = pnames)
    } else {
        if (length(notFounds))
            .warnnomatch("colnames", colmatches[notFounds])
        sampleMap[na.omit(colmatches), c("primary", "colname")]
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get the dimension names for
#' an \code{ExperimentList} using \code{\linkS4class{CharacterList}}
setMethod("dimnames", "ExperimentList", function(x) {
    list(IRanges::CharacterList(lapply(x, function(g) dimnames(g)[[1]])),
    IRanges::CharacterList(lapply(x, function(g) dimnames(g)[[2]])))
})

#' @describeIn MultiAssayExperiment Get the dimension names
#' for a \code{MultiAssayExperiment} object
setMethod("dimnames", "MultiAssayExperiment", function(x) {
    dimnames(experiments(x))
})

#' @export
.DollarNames.MultiAssayExperiment <- function(x, pattern = "")
    grep(pattern, names(colData(x)), value = TRUE)

#' @aliases $,MultiAssayExperiment-method
#' @exportMethod $
#' @rdname MultiAssayExperiment-methods
setMethod("$", "MultiAssayExperiment", function(x, name) {
    colData(x)[[name]]
})


.splitArgs <- function(args) {
    assayArgNames <- c("mcolname", "background", "type",
                       "make.names", "ranges")
    assayArgs <- args[assayArgNames]
    altArgs <- args[!names(args) %in% assayArgNames]
    assayArgs <- Filter(function(x) !is.null(x), assayArgs)
    list(assayArgs, altArgs)
}

.mergeMAE <- function(x, y) {
    if (any(names(x) %in% names(y)))
        stop("Provide unique experiment names")
    expz <- c(experiments(x), experiments(y))
    sampz <- rbind(sampleMap(x), sampleMap(y))
    coldx <- colData(x)
    coldx[["rnames"]] <- rownames(coldx)
    coldy <- colData(y)
    coldy[["rnames"]] <- rownames(coldy)
    cdatz <- S4Vectors::merge(coldx, coldy,
        by = intersect(names(coldx), names(coldy)), all = TRUE, sort = FALSE)
    rownames(cdatz) <- cdatz[["rnames"]]
    cdatz <- cdatz[, -which(names(cdatz) == "rnames")]
    metaz <- c(metadata(x), metadata(y))
    new("MultiAssayExperiment",
        ExperimentList = expz,
        colData = cdatz,
        sampleMap = sampz,
        metadata = metaz)
}

#' @describeIn MultiAssayExperiment Add a supported data class to the
#' \code{ExperimentList}
#'
#' @param sampleMap \code{c} method: a \code{sampleMap} \code{list} or
#' \code{DataFrame} to guide merge
#' @param mapFrom Either a \code{logical}, \code{character}, or \code{integer}
#' vector indicating the experiment(s) that have an identical colname order as
#' the experiment input(s)
#'
#' @examples
#' example("MultiAssayExperiment")
#'
#' ## Add an experiment
#' test1 <- mae[[1L]]
#' colnames(test1) <- rownames(colData(mae))
#'
#' ## Combine current MultiAssayExperiment with additional experiment
#' ## (no sampleMap)
#' c(mae, newExperiment = test1)
#'
#' test2 <- mae[[1L]]
#' c(mae, newExp = test2, mapFrom = 3L)
#'
setMethod("c", "MultiAssayExperiment",
    function(x, ..., sampleMap = NULL, mapFrom = NULL) {
    args <- list(...)
    if (!length(args))
        stop("Provide experiments or a 'MultiAssayExperiment' to concatenate")
    if (is(args[[1L]], "MultiAssayExperiment") && length(args) == 1L)
        return(.mergeMAE(x, args[[1L]]))
    exps <- ExperimentList(args)
    xmap <- sampleMap(x)
    cdata <- colData(x)
    if (!isEmpty(exps)) {
        if (!is.null(mapFrom)) {
            warning("Assuming column order in the data provided ",
                "\n matches the order in 'mapFrom' experiment(s) colnames")
            addMaps <- mapToList(sampleMap(x))[mapFrom]
            names(addMaps) <- names(exps)
            sampleMap <- mapply(function(x, y) {
                x[["colname"]] <- colnames(y)
                x
            }, addMaps, exps)
        } else if (is.null(sampleMap)) {
            sampleMap <- lapply(colnames(exps), .sampleMapFromExisting,
                colData = cdata, sampleMap = xmap)
            sampleMap <- listToMap(sampleMap)
        }
        if (is(sampleMap, "DataFrame") || is.data.frame(sampleMap))
            sampleMap <- mapToList(sampleMap)
        else if (!is.list(sampleMap))
            stop("'sampleMap' must be a 'DataFrame', 'data.frame', or 'list'")
        newListMap <- c(mapToList(xmap),
                        IRanges::SplitDataFrameList(sampleMap))
        x <- BiocGenerics:::replaceSlots(x,
            ExperimentList = c(experiments(x), exps),
            sampleMap = listToMap(newListMap)
        )
    }
    return(x)
})
