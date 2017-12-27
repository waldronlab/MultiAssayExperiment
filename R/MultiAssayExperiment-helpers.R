#' @include MultiAssayExperiment-methods.R
NULL

#' @name MultiAssayExperiment-helpers
#' @title A group of helper functions for manipulating and cleaning a
#' MultiAssayExperiment
#' @aliases intersectRows intersectColumns mergeReplicates replicated
#' complete.cases,MultiAssayExperiment-method
#' @description A set of helper functions were created to help clean and
#' manipulate a MultiAssayExperiment object. \code{intersectRows} also works
#' for \code{ExperimentList} objects.
#'
#' \itemize{
#'     \item complete.cases: Returns a logical vector corresponding to 'colData'
#'     rows that have data across all experiments
#'     \item isEmpty: Returns a logical \code{TRUE} value for zero length
#'     \code{MultiAssayExperiment} objects
#'     \item intersectRows: Takes all common rows across experiments,
#'     excludes experiments with empty rownames
#'     \item intersectColumns: A wrapper for \code{complete.cases} to return a
#'     \code{MultiAssayExperiment} with only those biological units that have
#'     measurements across all experiments
#'     \item replicated: A function that identifies multiple samples that
#'     originate from a single biological unit within each assay
#'     \item mergeReplicates: A function that combines replicated / repeated
#'     measurements across all experiments and is guided by the replicated
#'     return value
#'     \item longFormat: A \code{MultiAssayExperiment} method that
#'     returns a small and skinny \link{DataFrame}. The \code{colDataCols}
#'     arguments allows the user to append \code{colData} columns to the data.
#'     \item wideFormat: A function to return a wide \link{DataFrame} where
#'     each row represents an observation. Optional \code{colDataCols} can be
#'     added when using a \code{MultiAssayExperiment}.
#'     \item hasRowRanges: A function that identifies ExperimentList elements
#'     that have a \link[SummarizedExperiment]{rowRanges} method
#'     \item duplicated: (Deprecated) Returns a 'list' of 'LogicalList's that
#'     indicate what measurements originate from the same biological unit
#' }
#'
#' @param x A MultiAssayExperiment or ExperimentList
#'
#' @exportMethod complete.cases
setMethod("complete.cases", "MultiAssayExperiment", function(...) {
    args <- list(...)
    if (length(args) == 1L) {
        oldMap <- sampleMap(args[[1L]])
        listMap <- mapToList(oldMap)
        allPrimary <- Reduce(intersect,
            lapply(listMap, function(element) { element[["primary"]] }))
        rownames(colData(args[[1L]])) %in% allPrimary
    } else { stop("Provide only a 'MultiAssayExperiment'") }
})

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod isEmpty
setMethod("isEmpty", "MultiAssayExperiment", function(x) length(x) == 0L)

#' @rdname MultiAssayExperiment-helpers
#' @export
intersectRows <- function(x) {
    rows <- rownames(x)
    validRows <- Filter(length, rows)
    intRows <- Reduce(intersect, validRows)
    if (is(x, "MultiAssayExperiment"))
        x[intRows, , drop = FALSE]
    else if (is(x, "ExperimentList"))
        x[CharacterList(rep(list(intRows), length(validRows)))]
    else
        stop("Provide a valid class: ", class(x))
}

#' @rdname MultiAssayExperiment-helpers
#' @export
intersectColumns <- function(x) {
    comps <- complete.cases(x)
    x[, comps, drop = FALSE]
}

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("replicated", function(x) standardGeneric("replicated"))

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{replicated} function finds replicate samples in each
#' assay and returns a list of \linkS4class{LogicalList}s for each assay.
#' A \linkS4class{LogicalList} is given for each assay where each element
#' in such list corresponds to a biological unit entry in the \code{sampleMap},
#' (labeled as "primary"). Each element in this list of "primary" vectors is
#' a logical vector that identifies samples for that given biological unit.
#' Each logical vector is labeled with the "primary" identifier found in the
#' \code{sampleMap}. The \code{anyReplicated} will test to see if any of these
#' logical vectors has a summation value greater than one
#' (i.e., \strong{\code{sum(...) > 1L}}). These methods are not available for
#' an \code{ExperimentList} due to the unavailability of the \code{sampleMap}
#' structure.
setMethod("replicated", "MultiAssayExperiment", function(x) {
    listMap <- mapToList(sampleMap(x))
    lapply(listMap, function(assayDF) {
        pnames <- unique(assayDF[["primary"]])
        lmat <- vapply(pnames, function(x) {
            tots <- assayDF[["primary"]] %in% x
            if (sum(tots) <= 1L)
                tots <- rep(FALSE, nrow(assayDF))
            tots
        }, logical(nrow(assayDF)))
        resChunk <- LogicalList(lapply(seq_len(ncol(lmat)),
            function(x) lmat[, x]))
        names(resChunk) <- colnames(lmat)
        resChunk
    })
})

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("anyReplicated", function(x) standardGeneric("anyReplicated"))

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod anyReplicated
setMethod("anyReplicated", "MultiAssayExperiment", function(x) {
    reps <- replicated(x)
    vapply(reps, function(x) any(as.matrix(x)), logical(1L))
})

# mergeReplicates function ------------------------------------------------

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("mergeReplicates", function(x, replicates = list(),
                                       simplify = BiocGenerics::mean, ...)
    standardGeneric("mergeReplicates"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The \code{mergeReplicates} function is a house-keeping method
#' for a \code{MultiAssayExperiment} where only \code{complete.cases} are
#' returned, replicate measurements are averaged (by default), and columns are
#' aligned by the row order in \code{colData}. Additional arguments can be
#' passed on to the \code{simplify} function.
#'
#' @section mergeReplicates:
#' The \code{mergeReplicates} function makes use of the output from
#' \code{replicated} which will point out the duplicate measurements by
#' biological unit in the \code{MultiAssayExperiment}. This function will return
#' a \code{MultiAssayExperiment} with merged replicates. Additional arguments
#' can be provided to the simplify argument via the ellipsis (\ldots).
#'
#' @param replicates A list of \linkS4class{LogicalList}s
#' indicating multiple / duplicate entries for each biological unit per assay,
#' see \code{replicated} (default \code{replicated(x)}).
#' @param simplify A function for merging repeat measurements in experiments
#' as indicated by the \code{replicated} method for \code{MultiAssayExperiment}
#'
#' @exportMethod mergeReplicates
setMethod("mergeReplicates", "MultiAssayExperiment",
    function(x, replicates = replicated(x), simplify = BiocGenerics::mean, ...)
{
    if (!length(replicates))
        stop("'replicates' must be a list of technical replicates for each",
            "\n  biological unit. See '?replicated'.")
    experimentList <- mergeReplicates(
        x = experiments(x),
        replicates = replicates,
        simplify = simplify, ...)
    experiments(x) <- experimentList
    x
})

#' @describeIn ExperimentList Apply the mergeReplicates method on the
#' ExperimentList elements
#' @param replicates mergeReplicates: A \code{list} or \linkS4class{LogicalList}
#' where each element represents a sample and a vector of repeated measurements
#' for the sample
#' @param simplify A function for merging columns where duplicates are indicated
#' by replicates
#' @param ... Additional arguments. See details for more information.
setMethod("mergeReplicates", "ExperimentList",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...) {
        if (!length(replicates))
            stop("'replicates' must be a 'list' of replicated column elements",
                 "\n per biological unit")
        idx <- seq_along(x)
        names(idx) <- names(x)
        redList <- lapply(idx, function(i, element, simply,
                                        replicate, ...) {
            mergeReplicates(x = element[[i]], simplify = simply,
                            replicates = replicate[[i]], ...)
        }, element = x, simply = simplify,
        replicate = replicates, ...)
        ExperimentList(redList)
    })

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{mergeReplicates} "ANY" method consolidates duplicate
#' measurements for rectangular data structures, returns object of the same
#' class (endomorphic). The ellipsis or \code{\ldots} argument allows the
#' user to provide additional arguments to the \code{simplify} functional
#' argument.
setMethod("mergeReplicates", "ANY",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...) {
        object <- x
        if (is.list(replicates))
            replicates <- IRanges::LogicalList(replicates)
        if (is(object, "SummarizedExperiment"))
            x <- assay(object)
        if (is(object, "ExpressionSet"))
            x <- Biobase::exprs(object)
        if (any(any(replicates))) {
            uniqueCols <- apply(as.matrix(replicates), 2, function(cols) {
                !any(cols)
            })
            repeatList <- lapply(replicates, function(reps, rectangle) {
                if (any(reps)) {
                    repNames <- colnames(rectangle)[reps]
                    baseList <- as(split(rectangle[, repNames],
                                         seq_len(nrow(x))), "List")
                    result <- simplify(baseList, ...)
                    result <- matrix(result, ncol = 1,
                                     dimnames = list(NULL, repNames[[1L]]))
                    return(result)
                }
            }, rectangle = x)
            uniqueRectangle <- do.call(cbind, unname(repeatList))
            result <- cbind(uniqueRectangle, x[, uniqueCols, drop = FALSE])
            if (is(object, "SummarizedExperiment")) {
                # Keep only first replicate row in colData
                colDatIdx <- c(unname(min(which(replicates))),
                    which(uniqueCols))
                newColDat <- colData(object)[colDatIdx, , drop = FALSE]
                object <- initialize(object,
                    assays = Assays(SimpleList(result)), colData = newColDat)
            } else if (is(object, "ExpressionSet")) {
                # phenoData of ExpressionSet is lost
                object <- initialize(object, exprs = result)
            } else
                return(result)
        }
        return(object)
})

# longFormat function -----------------------------------------------------

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases longFormat
#' @section longFormat:
#' The longFormat method takes data from the \code{\link{ExperimentList}}
#' in a \code{\link{MultiAssayExperiment}} and returns a uniform
#' \code{\link{DataFrame}}. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' colData columns with the \code{colDataCols} argument
#' (\code{MultiAssayExperiment} method only). The \code{\ldots} argument
#' allows the user to specify the assay value for the
#' \linkS4class{SummarizedExperiment} assay function's \code{i} argument.
#'
#' @param object Any supported class object
#'
#' @export longFormat
setGeneric("longFormat", function(object, ...) standardGeneric("longFormat"))

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{longFormat} "ANY" class method, works with classes such as
#' \link{ExpressionSet} and \link{SummarizedExperiment} as well as \code{matrix}
#' to provide a consistent long and skinny \link{DataFrame}.
setMethod("longFormat", "ANY", function(object, ...) {
    rowNAMES <- rownames(object)
    nullROWS <- is.null(rowNAMES)
    args <- list(...)
    if (nullROWS)
        rowNAMES <- as.character(seq_len(nrow(object)))
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "SummarizedExperiment")) {
        object <- assay(object,
            i = if (!is.null(args[["i"]])) { args[["i"]] } else { 1L })
    }
    if (is(object, "matrix")) {
        object <- as.data.frame(object)
    }
    ## use stats::reshape instead of reshape2::melt
    if (nullROWS)
        rownames(object) <- rowNAMES
    object <- stats::reshape(object, idvar = "rowname",
        ids = rownames(object), times = names(object),
        timevar = "colname", varying = list(names(object)),
        direction = "long", v.names = "value")
    ## Reshape leaves rownames even if new.row.names = NULL
    rownames(object) <- NULL
    object <- object[, c("rowname", "colname", "value")]
    rectangle <- as(object, "DataFrame")
    rectangle[, "colname"] <- S4Vectors::Rle(rectangle[["colname"]])
    rectangle
})

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod longFormat
setMethod("longFormat", "ExperimentList", function(object, ...) {
    dataList <- lapply(seq_along(object), function(i, flatBox) {
        S4Vectors::DataFrame(assay = S4Vectors::Rle(names(object)[i]),
                             longFormat(flatBox[[i]], ...))
    }, flatBox = object)
    do.call(rbind, dataList)
})

#' @rdname MultiAssayExperiment-helpers
#' @param colDataCols selected colData columns to include in the output
setMethod("longFormat", "MultiAssayExperiment",
          function(object, colDataCols = NULL, ...) {
    addCols <- !is.null(colDataCols)
    longDataFrame <- longFormat(experiments(object), ...)
    primary <- S4Vectors::Rle(
        sampleMap(object)[match(longDataFrame[["colname"]],
                                sampleMap(object)[["colname"]]),
                          "primary"])
    longDataFrame <- S4Vectors::DataFrame(longDataFrame, primary = primary)
    longDataFrame <-
        longDataFrame[, c("assay", "primary", "rowname", "colname", "value")]
    if (addCols) {
        extraColumns <- colData(object)[, colDataCols, drop = FALSE]
        rowNameValues <- rownames(extraColumns)
        rownames(extraColumns) <- NULL
        matchIdx <- BiocGenerics::match(longDataFrame[["primary"]],
                                        rowNameValues)
        longDataFrame <- BiocGenerics::cbind(longDataFrame,
                                               extraColumns[matchIdx, ,
                                                            drop = FALSE])
    }
    return(longDataFrame)
})

# wideformat function -----------------------------------------------------

.uniteDF <- function(dframe, coi) {
    data.frame(
        feature = apply(dframe[, coi], 1L, paste, collapse = "_"),
        primary = dframe[["primary"]],
        value = dframe[["value"]]
    )
}

#' @rdname MultiAssayExperiment-helpers
#' @aliases wideFormat
#' @export
setGeneric("wideFormat", function(object, ...) standardGeneric("wideFormat"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @section wideFormat:
#' The \code{wideFormat} \code{MultiAssayExperiment} method returns standardized
#' wide \link{DataFrame} where each row represents an observation or biological
#' unit as represented in \code{colData}. Optionally, \code{colData} columns
#' can be added to the data output. The \code{wideFormat} method for an
#' \code{ExperimentList} returns a list of wideFormat \code{DataFrames}. The \code{\ldots}
#' argument allows the user to specify the assay number for
#' \linkS4class{SummarizedExperiment} assay extractor (i.e., \code{i} argument).
#' Additionally, the user may also specify \code{check.names} argument for the
#' resulting \code{DataFrame}. \strong{Note}. The "ANY" method returns a
#' slightly different wide format \code{DataFrame} due to missing biological
#' units information.
#'
#' @param key name of column whose values will used as variables in
#' the wide dataset, see the \code{timevar} argument in
#' \link[stats]{reshape}. If none are specified, assay, rowname, and colname
#' will be combined and named as the "feature" column (default "feature")
#' @param ... Additional arguments. See details.
#'
#' @exportMethod wideFormat
setMethod("wideFormat", "MultiAssayExperiment",
    function(object, colDataCols = NULL, key = "feature", ...)
{
    cnames <- colnames(object)
    check.names <- list(...)[["check.names"]]
    if (is.null(check.names)) check.names <- TRUE

    longDataFrame <- longFormat(object, colDataCols = colDataCols, ...)
    longDataFrame <- as.data.frame(longDataFrame)
    colsofinterest <- c("assay", "rowname")

    if (any(anyReplicated(object))) {
        dups <- replicated(object)
        lVects <- lapply(seq_along(dups), function(i, duplic) {
            assayname <- names(duplic)[[i]]
            logilist <- duplic[[i]]
            lmat <- as.matrix(logilist)
            rownames(lmat) <- names(logilist)
            colnames(lmat) <- cnames[[i]]
            lData <- longDataFrame[
                longDataFrame==assayname, c("primary", "colname")]
            apply(lData, 1L, function(x) lmat[x[1L], x[2L]])
        }, duplic = dups)

        longDataFrame[["replicated"]] <- unlist(lVects)
        splitDF <- split(longDataFrame, longDataFrame[["replicated"]])
        splitDF <- lapply(splitDF, function(x) x[, names(x) != "replicated"])

        splitDF <- lapply(names(splitDF), function(splitter) {
            dfchunk <- splitDF[[splitter]]
            colsofinterest <- if (splitter == "TRUE") {
                append(colsofinterest, "colname") } else { colsofinterest }
            .uniteDF(dfchunk, colsofinterest)
        })
        wideData <- do.call(rbind, splitDF)

    } else {
        wideData <- .uniteDF(longDataFrame, colsofinterest)
    }
    wideData <- stats::reshape(wideData, direction = "wide",
        idvar = "primary", timevar = key)
    names(wideData) <- gsub("value\\.", "", names(wideData))
    wideData <- wideData[
        match(rownames(colData(object)), wideData[["primary"]]), ]
    S4Vectors::DataFrame(wideData, check.names = check.names)
})

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The \code{wideFormat} method provides a representation of
#' \code{MultiAssayExperiment} where each row represents a particular
#' biological unit.
setMethod("wideFormat", "ANY", function(object, ...) {
    check.names <- list(...)[["check.names"]]
    if (is.null(check.names)) check.names <- TRUE
    args <- list(...)
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is.matrix(object) && is.null(rownames(object)))
        rownames(object) <- as.character(seq_len(object))
    if (!is.null(rownames(object)) && !is(object, "SummarizedExperiment"))
        object <- data.frame(rowname = rownames(object), object,
            stringsAsFactors = FALSE, check.names = check.names,
            row.names = NULL)
    if (is(object, "SummarizedExperiment")) {
        ## Ensure that rowData DataFrame has a rowname column
        ## Otherwise, use the rownames or first column
        rowDatNames <- names(rowData(object))
        rownameIn <- "rowname" %in% rowDatNames
        rowNAMES <- rownames(object)
        if (any(rownameIn)) {
            rowData(object) <- rowData(object)[rownameIn]
        } else if (!is.null(rowNAMES) || !length(rowDatNames)) {
            rowData(object) <- S4Vectors::DataFrame(rowname = rowNAMES)
        } else {
            warning("'rowname' column not in 'rowData' taking first one")
            rowData(object) <- rowData(object)[1L]
            names(rowData(object)) <- "rowname"
        }
        assayDat <- assay(object,
            i = if (!is.null(args[["i"]])) { args[["i"]] } else { 1L })
        object <- data.frame(rowname = rowData(object), assayDat,
            stringsAsFactors = FALSE, check.names = check.names,
            row.names = NULL)
    }
    S4Vectors::DataFrame(object)
})

.tryRowRanges <- function(obj) {
    res <- try(rowRanges(obj), silent = TRUE)
    if (!is(res, "try-error"))
        is(res, "GRanges")
    else
        FALSE
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases hasRowRanges
#' @section hasRowRanges:
#' The \code{hasRowRanges} method identifies assays with associated ranged
#' row data by directly testing the method on the object. The result from the
#' test must be a \linkS4class{GRanges} class object to satisfy the test.
#'
#' @export hasRowRanges
setGeneric("hasRowRanges", function(x) standardGeneric("hasRowRanges"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The \code{hasRowRanges} method identifies assays that support
#' a \link[SummarizedExperiment]{rowRanges} method \emph{and} return a
#' \linkS4class{GRanges} object.
setMethod("hasRowRanges", "MultiAssayExperiment", function(x) {
    hasRowRanges(experiments(x))
})

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod hasRowRanges
setMethod("hasRowRanges", "ExperimentList", function(x) {
    vapply(x, .tryRowRanges, logical(1L))
})

#' @rdname MultiAssayExperiment-helpers
#'
#' @param incomparables unused argument
#' @exportMethod duplicated
#' @aliases duplicated
#'
#' @details \strong{Deprecated:} For the \code{anyDuplicated} and
#' \code{duplicated} functions, the \code{incomparables} and ellipsis
#' \code{\ldots} arguments are not used. Neither \code{duplicated} nor
#' \code{anyDuplicated} is supported for \code{ExperimentList} due to an
#' unavailable \code{sampleMap}.
setMethod("duplicated", "MultiAssayExperiment",
          function(x, incomparables = FALSE, ...) {
    .Deprecated("replicated")
    replicated(x)
})
