#' @include MultiAssayExperiment-methods.R
NULL

#' @name MultiAssayExperiment-helpers
#' @title A group of helper functions for manipulating and cleaning a
#' MultiAssayExperiment
#' @aliases intersectRows intersectColumns mergeReplicates duplicated
#' @description A set of helper functions were created to help clean and
#' manipulate a MultiAssayExperiment object.
#'
#' \itemize{
#'     \item complete.cases: Returns a logical vector corresponding to 'colData'
#'     rows that have data across all experiments
#'     \item duplicated: Returns a 'list' of 'LogicalList's that indicate
#'     what measurements originate from the same biological unit
#'     \item intersectRows: Takes all common rows across experiments,
#'     excludes experiments with empty rownames
#'     \item intersectColumns: A wrapper for
#'     \link[=complete.cases,MultiAssayExperiment-method]{complete.cases} to
#'     return a MultiAssayExperiment with only those biological units that have
#'     measurements across all experiments
#'     \item mergeReplicates: A function that combines duplicated / repeated
#'     measurements across all experiments and is guided by the duplicated
#'     return value
#'     \item longFormat: A \code{MultiAssayExperiment} method that
#'     returns a small and skinny \link{DataFrame}. The \code{colDataCols}
#'     arguments allows the user to append \code{colData} columns to the data.
#'     \item wideFormat: A function to return a wide \link{DataFrame} where
#'     each row represents an observation. Optional \code{colDataCols} can be
#'     added when using a \code{MultiAssayExperiment}.
#' }
#'
#' @export
intersectRows <- function(x) {
    rows <- rownames(x)
    validRows <- Filter(length, rows)
    intRows <- Reduce(intersect, validRows)
    x[intRows, , drop = FALSE]
}

#' @rdname MultiAssayExperiment-helpers
#' @export
intersectColumns <- function(x) {
    comps <- complete.cases(x)
    x[, comps, drop = FALSE]
}

#' @rdname MultiAssayExperiment-helpers
#' @param x A MultiAssayExperiment
#' @param incomparables unused argument
#' @exportMethod duplicated
setMethod("duplicated", "MultiAssayExperiment",
          function(x, incomparables = FALSE, ...) {
    listMap <- mapToList(sampleMap(x))
    repList <- lapply(listMap, function(assayDF) {
        repeats <- unique(assayDF[["primary"]][
            duplicated(assayDF[["primary"]])])
        repSamps <- lapply(repeats, function(primary) {
            assayDF[["primary"]] %in% primary
        })
        names(repSamps) <- repeats
        repSamps
    })
    lapply(repList, IRanges::LogicalList)
})

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("mergeReplicates", function(x, replicates = list(),
                                       simplify = BiocGenerics::mean, ...)
    standardGeneric("mergeReplicates"))

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{mergeReplicates} function is a house-keeping method
#' for a \code{MultiAssayExperiment} where only \code{complete.cases} are
#' returned, replicate measurements are averaged (by default), and columns are
#' aligned by the row order in \code{colData}. Additional arguments can be
#' passed on to the \code{simplify} function.
#' @section mergeReplicates:
#' The \code{mergeReplicates} function makes use of the output from
#' \code{duplicated} which will point out the duplicate measurements by
#' biological unit in the \code{MultiAssayExperiment}. This function will return
#' a \code{MultiAssayExperiment} with merged replicates.
#' @param replicates A list of \linkS4class{LogicalList}s
#' indicating multiple / duplicate entries for each biological unit, see the
#' \code{duplicated} output
#' @param simplify A function for merging repeat measurements in experiments
#' as indicated by replicates for \code{MultiAssayExperiment}
#' @exportMethod mergeReplicates
setMethod("mergeReplicates", "MultiAssayExperiment",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...) {
        if (!length(replicates))
            replicates <- duplicated(x)
    experimentList <- mergeReplicates(x = experiments(x),
                                      replicates = replicates,
                                      simplify = simplify, ...)
    experiments(x) <- experimentList
    x
})

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{mergeReplicates} "ANY" method consolidates duplicate
#' measurements for rectangular data structures, returns object of the same
#' class (endomorphic)
setMethod("mergeReplicates", "ANY",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...) {
        object <- x
        if (is.list(replicates))
            replicates <- IRanges::LogicalList(replicates)
        if (is(object, "SummarizedExperiment"))
            x <- assay(object)
        if (is(object, "ExpressionSet"))
            x <- Biobase::exprs(object)
        if (length(replicates)) {
            uniqueCols <- apply(as.matrix(replicates), 2, function(cols) {
                !any(cols)
            })
            repeatList <- lapply(replicates, function(reps, rectangle) {
                if (length(reps)) {
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

#' @rdname MultiAssayExperiment-helpers
#' @aliases longFormat
#' @section longFormat:
#' The longFormat method takes data from the \code{\link{ExperimentList}}
#' in a \code{\link{MultiAssayExperiment}} and returns a uniform
#' \code{\link{DataFrame}}. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' colData columns with the \code{colDataCols} argument
#' (\code{MultiAssayExperiment} method only).
#'
#' @param object Any supported class object
#' @param ... Additional arguments. See details.
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
    if (nullROWS)
        rowNAMES <- rep(NA_character_, nrow(object))
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "SummarizedExperiment"))
        object <- assay(object)
    if (is(object, "matrix") && !nullROWS) {
        object <- reshape2::melt(object, varnames = c("rowname", "colname"),
                   as.is = TRUE)
    } else {
    object <- data.frame(rowname = rowNAMES, object,
                         stringsAsFactors = FALSE, check.names = FALSE,
                         row.names = NULL)
    object <- tidyr::gather(object, "colname", "value",
                            seq_along(object)[-1L])
    }
    rectangle <- S4Vectors::DataFrame(object)
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

#' @rdname MultiAssayExperiment-helpers
#' @aliases wideFormat
#' @export
setGeneric("wideFormat", function(object, ...) standardGeneric("wideFormat"))

#' @rdname MultiAssayExperiment-helpers
#' @section wideFormat:
#' The \code{wideFormat} \code{MultiAssayExperiment} method returns standardized
#' wide \link{DataFrame} where each row represents an observation or biological
#' unit as represented in \code{colData}. Optionally, \code{colData} columns
#' can be added to the data output. The \code{wideFormat} method for an
#' \code{ExperimentList} returns a list of wideFormat \code{DataFrames}. The
#' "ANY" method returns a wide format \code{DataFrame}.
#' @param key name of column whose values will used as variables in
#' the wide dataset from \link[tidyr]{spread}. If none are specified, assay,
#' rowname, and colname will be combined
setMethod("wideFormat", "MultiAssayExperiment",
    function(object, colDataCols = NULL, key = NULL, ...) {
        onetoone <- all(!lengths(duplicated(object)))
        longDataFrame <- longFormat(object, colDataCols = colDataCols, ...)
        longDataFrame <- as.data.frame(longDataFrame)
        if (is.null(key)) {
            if (onetoone) {
        longDataFrame <- tidyr::unite_(longDataFrame, "feature",
                                         c("assay", "rowname"))
        longDataFrame <- longDataFrame[, colnames(longDataFrame) != "colname"]
            } else {
        message("See ?mergeReplicates to combine replicated observations",
                "\n  to get one column per variable")
        longDataFrame <- tidyr::unite_(longDataFrame, "feature",
                                         c("assay", "rowname", "colname"))
            }
        wideDataFrame <- tidyr::spread_(longDataFrame, key = "feature",
                                         value = "value")
        } else {
        wideDataFrame <- tidyr::spread_(longDataFrame, key = key, value = "value")
        }
        wideDataFrame <- wideDataFrame[match(rownames(colData(object)),
            wideDataFrame[["primary"]]), ]
        wideDataFrame <- S4Vectors::DataFrame(wideDataFrame)
        return(wideDataFrame)
    })

#' @rdname MultiAssayExperiment-helpers
setMethod("wideFormat", "ExperimentList", function(object, ...) {
    lapply(object, wideFormat)
})

#' @rdname MultiAssayExperiment-helpers
setMethod("wideFormat", "ANY", function(object, ...) {
    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (!is.null(rownames(object)) && !is(object, "SummarizedExperiment"))
        object <- data.frame(rowname = rownames(object), object,
                             stringsAsFactors = FALSE, check.names = FALSE,
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
        object <- data.frame(rowname = rowData(object), object,
                             stringsAsFactors = FALSE, check.names = FALSE,
                             row.names = NULL)
    }
    S4Vectors::DataFrame(object)
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
            stop("'replicates' must be a 'list' of duplicated column elements",
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
