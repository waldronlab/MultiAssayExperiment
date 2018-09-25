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
#'     \item anyReplicated: Displays which assays have replicate measurements
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
#'     that have a \link[=RangedSummarizedExperiment-class]{rowRanges} method
#'     \item duplicated: (Deprecated) Returns a 'list' of 'LogicalList's that
#'     indicate what measurements originate from the same biological unit
#' }
#'
#' @param x A MultiAssayExperiment or ExperimentList
#' @param ... Additional arguments. See details for more information.
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

.matrixToList <- function(mat) {
    nout <- ncol(mat)
    listout <- vector("list", nout)
    for (i in seq_along(listout))
        listout[[i]] <- unname(mat[, i])
    names(listout) <- colnames(mat)
    listout
}

#' @rdname MultiAssayExperiment-helpers
#' @details The \code{replicated} function finds replicate measurements in each
#' assay and returns a list of \linkS4class{LogicalList}s.
#' Each element in a single \linkS4class{LogicalList} corresponds to a
#' biological or \emph{primary} unit as in the \code{sampleMap}. Below is a
#' small graphic for one particular biological unit in one assay, where the
#' logical vector corresponds to the number of measurements/samples in the
#' assay:
#' \preformatted{
#'  >      replicated(MultiAssayExperiment)
#'  (list str)       '-- $ AssayName
#'  (LogicalList str)      '-- [[ "Biological Unit" ]]
#'  Replicated if sum(...) > 1          '-- TRUE TRUE FALSE FALSE
#' }
#' \code{anyReplicated} determines if any of the assays have at least one
#' replicate. \emph{Note}. These methods are not available for the
#' \code{ExperimentList} class due to a missing \code{sampleMap} structure
#' (by design).
#'
#' @export
setMethod("replicated", "MultiAssayExperiment", function(x) {
    listMap <- mapToList(sampleMap(x))
    lapply(listMap, function(assayDF) {
        pnames <- unique(assayDF[["primary"]])
        lvect <- unlist(lapply(pnames, function(g) {
            tots <- assayDF[["primary"]] %in% g
            if (sum(tots) <= 1L)
                tots <- rep(FALSE, nrow(assayDF))
            tots
        }))
        lmat <- matrix(
            lvect, ncol = length(pnames), dimnames = list(NULL, pnames)
        )
        IRanges::LogicalList(.matrixToList(lmat))
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
setGeneric("mergeReplicates",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...)
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
    if (!length(replicates) || !identical(length(replicates), length(x)))
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
        expnames <- stats::setNames(names(x), names(x))
        redList <- lapply(expnames,
            function(i, element, simply, replicate, ...) {
                mergeReplicates(x = element[[i]], simplify = simply,
                    replicates = replicate[[i]], ...)
            }, element = x, simply = simplify, replicate = replicates, ...)
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
                colDatIdx <- c(unname(min(which(replicates[any(replicates)]))),
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

.longFormatANY <- function(object, i) {
    rowNAMES <- rownames(object)
    nullROWS <- is.null(rowNAMES)
    if (nullROWS)
        rowNAMES <- as.character(seq_len(nrow(object)))

    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "SummarizedExperiment"))
        object <- assay(object, i = i)
    if (is(object, "matrix"))
        object <- as.data.frame(object)

    ## use stats::reshape instead of reshape2::melt
    if (nullROWS)
        rownames(object) <- rowNAMES
    object <- stats::reshape(object, idvar = "rowname",
        ids = rownames(object), times = names(object),
        timevar = "colname", varying = list(names(object)),
        direction = "long", v.names = "value")
    ## Reshape leaves rownames even if new.row.names = NULL
    rownames(object) <- NULL
    object[, c("rowname", "colname", "value")]
}

.longFormatElist <- function(object, i) {
    if (!is(object, "ExperimentList"))
        stop("<internal> Not an 'ExperimentList' input")
    objnames <- stats::setNames(names(object), names(object))
    lapply(objnames, function(nameidx, flatBox) {
        data.frame(assay = nameidx,
            .longFormatANY(flatBox[[nameidx]], i = i),
            stringsAsFactors = FALSE)
        }, flatBox = object)
}

.matchAddColData <- function(reshaped, colData, colDataCols) {
    extraColumns <- as.data.frame(colData[, colDataCols, drop = FALSE])
    rowNameValues <- rownames(extraColumns)
    rownames(extraColumns) <- NULL
    matchIdx <- match(reshaped[["primary"]], rowNameValues)
    matchedDF <- extraColumns[matchIdx, , drop = FALSE]
    rownames(matchedDF) <- NULL
    cbind(reshaped, matchedDF)
}

.mapOrderPrimary <- function(flatbox, samplemap) {
    primary <- samplemap[match(flatbox[["colname"]], samplemap[["colname"]]),
        "primary"]

    reshaped <- data.frame(flatbox, primary = primary,
        stringsAsFactors = FALSE)
    reshaped[, c("assay", "primary", "rowname", "colname", "value")]
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases longFormat
#'
#' @details The \code{longFormat} "ANY" class method, works with classes such as
#' \link{ExpressionSet} and \link{SummarizedExperiment} as well as \code{matrix}
#' to provide a consistent long and skinny \link{DataFrame}.
#'
#' @section longFormat:
#' The longFormat method takes data from the \code{\link{ExperimentList}}
#' in a \code{\link{MultiAssayExperiment}} and returns a uniform
#' \code{\link{DataFrame}}. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' colData columns with the \code{colDataCols} argument
#' (\code{MultiAssayExperiment} method only). The \code{i} argument
#' allows the user to specify the assay value in a
#' \linkS4class{SummarizedExperiment}. It directly relates to the \code{i}
#' argument in the assay method.
#'
#' @param object Any supported class object
#' @param colDataCols A \code{character}, \code{logical}, or \code{numeric}
#' index for \code{colData} columns to be included
#' @param i The assay indicator for \linkS4class{SummarizedExperiment}
#' objects (default 1L)
#'
#' @export longFormat
longFormat <- function(object, colDataCols = NULL, i = 1L) {
    if (is(object, "ExperimentList"))
        return(do.call(rbind, .longFormatElist(object, i = i)))
    else if (!is(object, "MultiAssayExperiment"))
        return(.longFormatANY(object, i = i))

    longDataFrame <- do.call(function(...) rbind(..., make.row.names = FALSE),
        .longFormatElist(experiments(object), i = i))

    longDataFrame <- .mapOrderPrimary(longDataFrame, sampleMap(object))

    if (!is.null(colDataCols))
        longDataFrame <-
            .matchAddColData(longDataFrame, colData(object), colDataCols)

    as(longDataFrame, "DataFrame")
}

# wideformat function -----------------------------------------------------

.wideFormatANY <- function(object, i, check.names) {
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
        assayDat <- assay(object, i = i)
        object <- data.frame(rowname = rowData(object), assayDat,
            stringsAsFactors = FALSE, check.names = check.names,
            row.names = NULL)
    }
    object
}

.metadataCOLS <- function(metcols, collapser, coldatcols) {
    namesList <- lapply(strsplit(metcols, collapser), `[`)
    sqnames <- lapply(namesList, function(nam)
        append(nam, rep(NA_character_, 3L - length(nam))))
    doFrame <- do.call(rbind.data.frame, sqnames)
    names(doFrame) <- c("sourceName", "rowname", "colname")
    doFrame[["sourceName"]] <-
        gsub("primary", "colDataRows", doFrame[["sourceName"]])
    doFrame[doFrame[["sourceName"]] %in% coldatcols, "sourceName"] <-
        "colDataCols"
    doFrame
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @section wideFormat:
#' The \code{wideFormat} function returns standardized wide \link{DataFrame}
#' where each row represents a biological unit as in the \code{colData}.
#' Depending on the data and setup, biological units can be patients, tumors,
#' specimens, etc. Optionally, \code{colData} columns can be added to the
#' wide data output (see the \code{colDataCols} argument). Metadata columns are
#' generated based on the names produced in the wide format
#' \linkS4class{DataFrame}. These can be accessed via the
#' \link[=Vector-class]{mcols} function. See the \code{Arguments} and
#' \code{longFormat} sections for argument descriptions.
#'
#' @param check.names (logical default TRUE) Column names of the output
#' \code{DataFrame} will be checked for syntactic validity and made unique,
#' if necessary
#' @param collapse (character default "_") A single string delimiter for output
#' column names. In \code{wideFormat}, experiments and rownames (and when
#' replicate samples are present, colnames) are seperated by this delimiter
#'
#' @export wideFormat
wideFormat <- function(object, colDataCols = NULL, check.names = TRUE,
    collapse = "_", i = 1L) {

    collSymbol <- "///"
    key <- "feature"

    if (is.null(colDataCols)) colDataCols <- character(0L)
    nameFUN <- if (check.names) make.names else I
    cnames <- colnames(object)
    longList <- .longFormatElist(experiments(object), i = i)
    longList <- lapply(longList, .mapOrderPrimary, sampleMap(object))
    colsofinterest <- c("assay", "rowname")

    anyReps <- anyReplicated(object)
    if (any(anyReps)) {
        dups <- replicated(object)[anyReps]
        indx <- names(which(anyReps))
        lVects <- lapply(indx, function(expname, duplic) {
            logilist <- duplic[[expname]]
            lmat <- as.matrix(logilist)
            rownames(lmat) <- names(logilist)
            colnames(lmat) <- cnames[[expname]]
            lData <- longList[[expname]][, c("primary", "colname")]
            apply(lData, 1L, function(x) lmat[x[1L], x[2L]])
        }, duplic = dups)

        repList <- Map(function(x, y) { x[y, , drop = FALSE] },
            x = longList[anyReps], y = lVects)

        longList[anyReps] <- Map(function(x, y) { x[!y, , drop = FALSE] },
            x = longList[anyReps], y = lVects)

        longList <- lapply(longList, function(x)
            tidyr::unite(x[, names(x) != "colname"], col = !!key, colsofinterest,
                sep = collSymbol))

        repList <- lapply(repList, function(x)
            tidyr::unite(x, col = !!key, c(colsofinterest, "colname"),
                sep = collSymbol))

        wideData <- c(longList, repList)
    } else {

        wideData <- lapply(longList, function(x)
            tidyr::unite(x[, names(x) != "colname"], col = !!key,
                colsofinterest, sep = collSymbol))
    }

    wideData <- lapply(wideData, function(flox) {
        flox <- tidyr::spread(flox, key = {key}, value = "value")
        .matchAddColData(flox, colData(object), colDataCols)
    })
    wideDF <- Reduce(function(x, y)
        merge(x, y, by = intersect(names(x), names(y)), all = TRUE), wideData)
    wideDF <- as(wideDF, "DataFrame")

    metadat <- .metadataCOLS(names(wideDF), collSymbol, colDataCols)
    mcols(wideDF) <- metadat
    names(wideDF) <- nameFUN(gsub(collSymbol, collapse, names(wideDF)))

    wideDF[order(match(wideDF[["primary"]], rownames(colData(object)))), ]
}

# hasRowRanges section ----------------------------------------------------

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
#' a \link[=RangedSummarizedExperiment-class]{rowRanges} method \emph{and}
#' return a \linkS4class{GRanges} object.
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
