#' @include MultiAssayExperiment-subset.R
NULL

.warnnomatch <- function(type, values) {
    warning("'", type, "': ",
        paste(S4Vectors:::selectSome(values), collapse = ", "),
            "\n  could not be matched")
}

.generateMiniMap <- function(colData, sampleMap, colnames) {
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
#' test1 <- myMultiAssayExperiment[[1L]]
#' colnames(test1) <- rownames(colData(myMultiAssayExperiment))
#'
#' ## Combine current MultiAssayExperiment with additional experiment
#' ## (no sampleMap)
#' c(myMultiAssayExperiment, newExperiment = test1)
#'
#' test2 <- myMultiAssayExperiment[[1L]]
#' c(myMultiAssayExperiment, newExp = test2, mapFrom = 3L)
#'
setMethod("c", "MultiAssayExperiment",
    function(x, ..., sampleMap = NULL, mapFrom = NULL) {
    exps <- ExperimentList(...)
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
            sampleMap <- lapply(colnames(exps), .generateMiniMap,
                colData = cdata, sampleMap = xmap)
            sampleMap <- listToMap(sampleMap)
        }
        if (is(sampleMap, "DataFrame") || is.data.frame(sampleMap))
            sampleMap <- mapToList(sampleMap)
        else if (!is.list(sampleMap))
            stop("'sampleMap' must be a 'DataFrame', 'data.frame', or 'list'")
        newListMap <- c(mapToList(xmap),
                        IRanges::SplitDataFrameList(sampleMap))
        sampleMap(x) <- listToMap(newListMap)
        experiments(x) <- c(experiments(x), exps)
        validObject(x)
    }
    return(x)
})
