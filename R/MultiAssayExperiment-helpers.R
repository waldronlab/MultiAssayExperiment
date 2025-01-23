#' @include MultiAssayExperiment-methods.R
NULL

#' @name MultiAssayExperiment-helpers
#'
#' @title A group of helper functions for manipulating and cleaning a
#' MultiAssayExperiment
#'
#' @aliases intersectRows intersectColumns mergeReplicates replicated
#' complete.cases,MultiAssayExperiment-method
#'
#' @description A set of helper functions were created to help clean and
#' manipulate a MultiAssayExperiment object. `intersectRows` also works
#' for `ExperimentList` objects.
#'
#' * complete.cases: Returns a logical vector corresponding to 'colData'
#'     rows that have data across all experiments
#' * isEmpty: Returns a logical `TRUE` value for zero length
#'     `MultiAssayExperiment` objects
#' * intersectRows: Takes all common rows across experiments,
#'     excludes experiments with empty rownames
#' * intersectColumns: A wrapper for `complete.cases` to return a
#'     `MultiAssayExperiment` with only those biological units that have
#'     measurements across all experiments
#' * replicated: Identifies, via logical vectors, `colname`s that
#'     originate from a single biological unit within each assay
#' * replicates: Provides the replicate `colname`s found with
#'     the `replicated` function by their name, empty list if none
#' * anyReplicated: Whether the assay has replicate measurements
#' * showReplicated: Displays the actual columns that are replicated per
#'     assay and biological unit, i.e., `primary` value (`colData`
#'     rowname) in the `sampleMap`
#' * mergeReplicates: A function that combines replicated / repeated
#'     measurements across all experiments and is guided by the replicated
#'     return value
#' * longForm: A `MultiAssayExperiment` method that returns a small and
#'     skinny [`DataFrame`][S4Vectors::DataFrame-class]. The `colDataCols`
#'     arguments allows the user to append `colData` columns to the data.
#' * wideFormat: A function to reshape the data in a
#'     `MultiAssayExperiment` to a "wide" format
#'     [`DataFrame`][S4Vectors::DataFrame-class]. Each row in the `DataFrame`
#'     represents an observation (corresponding to an entry in the `colData`).
#'     If replicates are present, their data will be appended at
#'     the end of the corresponding row and will generate additional `NA` data.
#'     It is recommended to remove or consolidate technical replicates with
#'     `mergeReplicates`. Optional `colDataCols` can be added when the
#'     original object is a `MultiAssayExperiment`.
#' * hasRowRanges: A function that identifies ExperimentList elements
#'     that have a
#'     [`rowRanges`][SummarizedExperiment::RangedSummarizedExperiment-class]
#'     method
#' * hasRowData: A function that identifies ExperimentList elements
#'     that have a
#'     [`rowData`][SummarizedExperiment::SummarizedExperiment-class]
#'     method
#' * getWithColData: A convenience function for extracting an assay
#'     and associated colData
#' * renamePrimary: A convenience function to rename the primary
#'     biological units as represented in the `rownames(colData)`
#' * renameColname: A convenience function to rename the colnames
#'     of a particular assay
#'
#' @param x A MultiAssayExperiment or ExperimentList
#'
#' @param ... Additional arguments. See details for more information.
#'
#' @return See the itemized list in the description section for details.
#'
#' @examples
#'
#' example(MultiAssayExperiment)
#'
#' complete.cases(mae)
#'
#' isEmpty(MultiAssayExperiment())
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
#'
#' @details The `replicated` function finds replicate measurements in each
#' assay and returns a list of [`LogicalList`][IRanges::AtomicList]s.
#' Each element in a single [`LogicalList`][IRanges::AtomicList] corresponds to
#' a biological or _primary_ unit as in the `sampleMap`. Below is a
#' small graphic for one particular biological unit in one assay, where the
#' logical vector corresponds to the number of measurements/samples in the
#' assay:
#' \preformatted{
#'  >      replicated(MultiAssayExperiment)
#'  (list str)       '-- $ AssayName
#'  (LogicalList str)      '-- [[ "Biological Unit" ]]
#'  Replicated if sum(...) > 1          '-- TRUE TRUE FALSE FALSE
#' }
#' `anyReplicated` determines if any of the assays have at least one
#' replicate. _Note_. These methods are not available for the
#' `ExperimentList` class due to a missing `sampleMap` structure
#' (by design).
#' `showReplicated` returns a list of [`CharacterList`][IRanges::AtomicList]s
#' where each element corresponds to the the biological or _primary_ units that
#' are replicated in that assay element. The values in the inner list are
#' the `colnames` in the assay that are technical replicates.
#'
#' @export
setMethod("replicated", "MultiAssayExperiment", function(x) {
    listMap <- mapToList(sampleMap(x))
    lapply(listMap, function(assayDF) {
        pnames <- assayDF[["primary"]]
        IRanges::LogicalList(lapply(
            S4Vectors::splitAsList(pnames, pnames),
            function(g) {
                if (identical(length(g), 1L))
                    rep(FALSE, length(pnames))
                else
                    pnames %in% g
            }
        ))
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

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("showReplicated", function(x) standardGeneric("showReplicated"))

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod showReplicated
setMethod("showReplicated", "MultiAssayExperiment", function(x) {
    clnames <- Filter(length, colnames(x))
    replicates <- replicated(x)[names(clnames)]
    Map(
        function(y, z) {
            IRanges::CharacterList(
                Filter(length, lapply(z, function(g) y[g]))
            )
        },
        y = clnames,
        z = replicates
    )
})

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("replicates", function(x, ...) standardGeneric("replicates"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The `replicates` function (noun) returns the `colname`s
#'   from the `sampleMap` that were identified as replicates. It returns a
#'   list of [`CharacterList`][IRanges::AtomicList]s for each assay present in
#'   the `MultiAssayExperiment` and an inner entry for each biological unit that
#'   has replicate observations in that assay.
#'
#' @export
setMethod("replicates", "MultiAssayExperiment", function(x, ...) {
    listMap <- mapToList(sampleMap(x))
    lapply(
        X = listMap,
        FUN = function(assayDF) {
            Filter(
                f = function(y) {
                    length(y) > 1
                },
                x = S4Vectors::splitAsList(
                    assayDF[["colname"]], assayDF[["primary"]]
                )
            )
        }
    )
})

# mergeReplicates function ------------------------------------------------

#' @rdname MultiAssayExperiment-helpers
#' @export
setGeneric("mergeReplicates",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...)
        standardGeneric("mergeReplicates"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The `mergeReplicates` function is a house-keeping method
#' for a `MultiAssayExperiment` where only `complete.cases` are
#' returned. This by-assay operation averages replicate measurements
#' (by default) and columns are aligned by the row order in `colData`.
#' Users can provide their own function for merging replicates with the
#' `simplify` functional argument. Additional inputs `...` are
#' sent to the 'simplify' function.
#'
#' @section mergeReplicates:
#' The `mergeReplicates` function makes use of the output from
#' `replicated` which will point out the duplicate measurements by
#' biological unit in the `MultiAssayExperiment`. This function will
#' return a `MultiAssayExperiment` with merged replicates. Additional
#' arguments can be provided to the simplify argument via the ellipsis
#' (`...`). For example, when replicates "TCGA-B" and "TCGA-A" are found in
#' the assay, the name of the first appearing replicate is taken (i.e., "B").
#' Note that a typical use case of merging replicates occurs when there are
#' multiple measurements on the **same** sample (within the same assay)
#' and can therefore be averaged.
#'
#' @param replicates A list of [`LogicalList`][IRanges::AtomicList]s
#' indicating multiple / duplicate entries for each biological unit per assay,
#' see `replicated` (default `replicated(x)`).
#' @param simplify A function for merging repeat measurements in experiments
#' as indicated by the `replicated` method for `MultiAssayExperiment`
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
#'
#' @param replicates mergeReplicates: A `list` or
#'   [`LogicalList`][IRanges::AtomicList] where each element represents a
#'   sample and a vector of repeated measurements for the sample
#'
#' @param simplify A function for merging columns where duplicates are indicated
#' by replicates
#'
#' @param ... Additional arguments. See details for more information.
setMethod("mergeReplicates", "ExperimentList",
    function(x, replicates = list(), simplify = BiocGenerics::mean, ...) {
        if (!length(replicates))
            stop("'replicates' must be a 'list' of replicated column elements",
                 "\n per biological unit")
        expnames <- .setNames(names(x), names(x))
        redList <- lapply(expnames,
            function(i, element, simply, replicate, ...) {
                mergeReplicates(x = element[[i]], simplify = simply,
                    replicates = replicate[[i]], ...)
            }, element = x, simply = simplify, replicate = replicates, ...)
        ExperimentList(redList)
    }
)

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The `mergeReplicates` "ANY" method consolidates duplicate
#' measurements for rectangular data structures, returns object of the same
#' class (endomorphic). The ellipsis or `...` argument allows the
#' user to provide additional arguments to the `simplify` functional
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

# longFrom method ---------------------------------------------------------

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The `longForm` "ANY" class method, works with classes such as
#' [`ExpressionSet`][Biobase::ExpressionSet] and
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class] as
#' well as `matrix` to provide a consistent long and skinny
#' [`DataFrame`][S4Vectors::DataFrame-class].
#'
#' @section longForm:
#' The 'longForm' method takes data from the [`ExperimentList`]
#' in a `MultiAssayExperiment` and returns a uniform
#' `DataFrame`. The resulting DataFrame has columns indicating
#' primary, rowname, colname and value. This method can optionally include
#' columns of the MultiAssayExperiment colData named by `colDataCols` character
#' vector argument. (`MultiAssayExperiment` method only). The `i` argument
#' allows the user to specify the assay value for the
#' `SummarizedExperiment` assay function's `i` argument.
#'
#' @param object Any supported class object
#'
#' @param colDataCols A `character`, `logical`, or `numeric`
#'     index for `colData` columns to be included
#'
#' @param i longForm: The i-th assay in
#'     `SummarizedExperiment`-like objects. A vector input is
#'     supported in the case that the `SummarizedExperiment` object(s) has more
#'     than one assay (default 1L),
#'     renameColname: Either a `numeric` or `character` index
#'     indicating the assay whose colnames are to be renamed
#'
#' @importFrom BiocBaseUtils lifeCycle
#' @importFrom BiocGenerics longForm
#'
#' @exportMethod longForm
setMethod(
    "longForm", "MultiAssayExperiment",
    function(object, colDataCols = NULL, i = 1L, ...) {
        if (any(.emptyAssays(experiments(object))))
            object <- .dropEmpty(object, warn = FALSE)

        longDataFrame <- longForm(experiments(object), i = i)

        longDataFrame <- .mapOrderPrimary(longDataFrame, sampleMap(object))

        if (!is.null(colDataCols))
            longDataFrame <-
                .matchAddColData(longDataFrame, colData(object), colDataCols)

        as(longDataFrame, "DataFrame")
    }
)

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod longForm
setMethod(
    "longForm", "ExperimentList",
    function(object, colDataCols, i = 1L, ...) {
        samelength <- identical(length(object), length(i))
        if (!samelength && identical(length(i), 1L))
            i <- rep(i, length(object))
        res <- mapply(
            function(obj, obname, idx) {
                data.frame(
                    assay = obname,
                    longForm(obj, i = idx),
                    stringsAsFactors = FALSE
                )
            }, obj = object, obname = names(object), idx = i, SIMPLIFY = FALSE
        )
        do.call(
            function(...) rbind(..., make.row.names = FALSE),
            res
        )
    }
)

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod longForm
setMethod("longForm", "ANY", function(object, colDataCols, i = 1L, ...) {
    rowNAMES <- rownames(object)
    if (is.null(rowNAMES)) rowNames <- as.character(seq_len(nrow(object)))

    if (is(object, "ExpressionSet"))
        object <- Biobase::exprs(object)
    if (is(object, "SummarizedExperiment") || is(object, "RaggedExperiment"))
        object <- assay(object, i = i)

    BiocBaseUtils::checkInstalled("reshape2")

    res <- reshape2::melt(
        object, varnames = c("rowname", "colname"), value.name = "value"
    )
    if (!is.character(res[["rowname"]]))
        res[["rowname"]] <- as.character(res[["rowname"]])
    res
})

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
#' The `wideFormat` function returns standardized wide [`DataFrame`]
#' where each row represents a biological unit as in the `colData`.
#' Depending on the data and setup, biological units can be patients, tumors,
#' specimens, etc. Metadata columns are
#' generated based on the names produced in the wide format
#' `DataFrame`. These can be accessed via the
#' [`mcols()`][S4Vectors::Vector-class] function.
#' See the `wideFormat` section for description of the `colDataCols` and
#' `i` arguments.
#'
#' @param check.names (logical default TRUE) Column names of the output
#' `DataFrame` will be checked for syntactic validity and made unique,
#' if necessary
#'
#' @param collapse (character default "_") A single string delimiter for output
#' column names. In `wideFormat`, experiments and rownames (and when
#' replicate samples are present, colnames) are seperated by this delimiter
#'
#' @export wideFormat
wideFormat <- function(object, colDataCols = NULL, check.names = TRUE,
    collapse = "_", i = 1L) {

    collSymbol <- "///"
    key <- "feature"

    if (any(.emptyAssays(experiments(object))))
        object <- .dropEmpty(object, warn = FALSE)

    if (is.null(colDataCols)) colDataCols <- character(0L)
    nameFUN <- if (check.names) make.names else I
    cnames <- colnames(object)
    longDataFrame <- longForm(experiments(object), i = i)
    longList <- split(longDataFrame, longDataFrame[["assay"]])
    longList <- lapply(longList, .mapOrderPrimary, sampleMap(object))
    colsofinterest <- c("assay", "rowname")

    anyReps <- anyReplicated(object)
    if (any(anyReps)) {
        indx <- names(which(anyReps))
        dups <- replicated(object)[indx]
        lVects <- lapply(indx, function(expname, duplic) {
            logilist <- duplic[[expname]]
            lmat <- as.matrix(logilist)
            rownames(lmat) <- names(logilist)
            colnames(lmat) <- cnames[[expname]]
            lData <- longList[[expname]][, c("primary", "colname")]
            apply(lData, 1L, function(x) lmat[x[1L], x[2L]])
        }, duplic = dups)

        repList <- Map(function(x, y) { x[y, , drop = FALSE] },
            x = longList[indx], y = lVects)

        longList[indx] <- Map(function(x, y) { x[!y, , drop = FALSE] },
            x = longList[indx], y = lVects)

        longList <- lapply(longList, function(x)
            tidyr::unite(x[, names(x) != "colname"], col = !!key,
                colsofinterest, sep = collSymbol)
        )

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
        flox <- tidyr::pivot_wider(flox, names_from = key)
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
#'
#' @section hasRowRanges:
#' The `hasRowRanges` method identifies assays with associated ranged
#' row data by directly testing the method on the object. The result from the
#' test must be a [`GRanges`][GenomicRanges::GRanges-class] class object to
#' satisfy the test.
#'
#' @export hasRowRanges
setGeneric("hasRowRanges", function(x) standardGeneric("hasRowRanges"))

#' @rdname MultiAssayExperiment-helpers
#'
#' @details The `hasRowRanges` method identifies assays that support
#' a [`rowRanges`][SummarizedExperiment::RangedSummarizedExperiment-class]
#' method _and_ return a [`GRanges`][GenomicRanges::GRanges-class] object.
setMethod("hasRowRanges", "MultiAssayExperiment", function(x) {
    hasRowRanges(experiments(x))
})

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod hasRowRanges
setMethod("hasRowRanges", "ExperimentList", function(x) {
    vapply(x, .tryRowRanges, logical(1L))
})

# hasRowData section ------------------------------------------------------

#' @importFrom SummarizedExperiment rowData
.tryRowData <- function(obj) {
    res <- try(rowData(obj), silent = TRUE)
    if (!is(res, "try-error"))
        is(res, "DataFrame")
    else
        FALSE
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases hasRowData
#'
#' @details The `hasRowData` method identifies assays that support a
#'   [`rowData`][SummarizedExperiment::SummarizedExperiment-class] method _and_
#'   return a [`DataFrame`][S4Vectors::DataFrame-class] object.
#'
#' @export
setGeneric("hasRowData", function(x) standardGeneric("hasRowData"))

#' @describeIn MultiAssayExperiment-helpers The `hasRowData` method identifies
#'   experiments that have a `rowData` method via direct testing
#'
#' @exportMethod hasRowData
setMethod("hasRowData", "MultiAssayExperiment", function(x) {
    hasRowData(experiments(x))
})

#' @describeIn MultiAssayExperiment-helpers The `hasRowData` method identifies
#'   experiments that have a `rowData` method via direct testing
#'
#' @exportMethod hasRowData
setMethod("hasRowData", "ExperimentList", function(x) {
    vapply(x, .tryRowData, logical(1L))
})

#' @rdname MultiAssayExperiment-helpers
#'
#' @param mode String indicating how `MultiAssayExperiment`
#'     column-level metadata should be added to the
#'     `SummarizedExperiment` `colData`.
#'
#' @param verbose `logical(1)` Whether to `suppressMessages` on subsetting
#'     operations in `getWithColData` (default FALSE)
#'
#' @aliases getWithColData
#'
#' @section getWithColData:
#' The `getWithColData` function allows the user to conveniently extract
#' a particular assay as indicated by the **`i`** index argument. It
#' will also attempt to provide the
#' [`colData`][SummarizedExperiment::SummarizedExperiment-class]
#' along with the extracted object using the `colData<-` replacement
#' method when possible. Typically, this method is available for
#' [`SummarizedExperiment`][SummarizedExperiment::SummarizedExperiment-class]
#' and `RaggedExperiment` classes.
#'
#' The setting of `mode` determines how the `colData`
#' is added. If `mode="append"`, the `MultiAssayExperiment`
#' metadata is appended onto that of the `SummarizedExperiment`.
#' If any fields are duplicated by name, the values in the
#' `SummarizedExperiment` are retained, with a warning emitted if
#' the values are different.  For `mode="replace"`, the
#' `MultiAssayExperiment` metadata replaces that of the
#' `SummarizedExperiment`, while for `mode="none"`,
#' no replacement or appending is performed.
#'
#' @export getWithColData
getWithColData <- function(x, i, mode=c("append", "replace"), verbose = FALSE) {
    if (!is(x, "MultiAssayExperiment"))
        stop("Provide a MultiAssayExperiment as input")

    stopifnot(is.numeric(i) || is.character(i),
        identical(length(i), 1L), !is.na(i), !is.logical(i),
        is.character(mode), !is.na(mode), !is.logical(mode))

    FUN <- if (!verbose) suppressMessages else force
    mae <- FUN(x[, , i, drop = FALSE])
    prims <- sampleMap(mae)[["primary"]]
    if (anyDuplicated(prims))
        warning("Duplicating colData rows due to replicates in 'replicated(x)'",
            call. = FALSE)
    expanded <- colData(mae)[prims, , drop = FALSE]
    exObj <- mae[[1L]]

    tryCatch({
        existing <- colData(exObj)
    }, error = function(e) {
        stop(
            "Extracted class does not support 'colData':",
            "\n    ", conditionMessage(e), call. = FALSE
        )
    })

    mode <- match.arg(mode)
    if (identical(mode, "replace") || !length(existing))
            colData(exObj) <- expanded
    else if (identical(mode, "append")) {
        # Prune out duplicate entries
        common <- intersect(colnames(existing), colnames(expanded))
        if (length(common)) {
            warnid <- vapply(common, function(varname, x, y) {
                !identical(x[[varname]], y[[varname]])
            }, logical(1L), x = existing, y = expanded)
            .warning(
                "Ignoring redundant column names in 'colData(x)': ",
                paste(common[warnid], collapse = ", ")
            )
        }
        leftovers <- expanded[, setdiff(colnames(expanded), common), drop=FALSE]
        colData(exObj) <- cbind(existing, leftovers)
    }

    exObj
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases renamePrimary
#'
#' @section rename*:
#' The `renamePrimary` function allows the user to conveniently change the
#' actual names of the primary biological units as seen in
#' `rownames(colData)`. `renameColname` allows the user to change the
#' names of a particular assay based on index `i`. `i` can either be
#' a single numeric or character value. See `colnames<-` method for
#' renaming multiple colnames in a `MultiAssayExperiment`.
#'
#' @param value renamePrimary: A `character` vector of the same length as
#' the existing `rownames(colData)` to use for replacement,
#' renameColname: A `CharacterList` or `list` with matching
#' `lengths` to replace `colnames(x)`
#'
#' @examples
#'
#' ## renaming biological units (primary)
#'
#' mae2 <- renamePrimary(mae, paste0("pt", 1:4))
#' colData(mae2)
#' sampleMap(mae2)
#'
#' @export renamePrimary
renamePrimary <- function(x, value) {
    coldata <- colData(x)
    old <- rownames(coldata)
    if (length(old) != length(value))
        stop("'value' length does not match 'length(rownames(colData))'")
    samplemap <- sampleMap(x)
    nprime <- value[match(samplemap[["primary"]], old)]
    samplemap[["primary"]] <- nprime
    rownames(coldata) <- value
    BiocBaseUtils::setSlots(
        object = x,
        sampleMap = samplemap,
        colData = coldata
    )
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases renameColname
#'
#' @examples
#'
#' ## renaming observational units (colname)
#'
#' mae2 <- renameColname(mae, i = "Affy", paste0("ARRAY", 1:4))
#' colnames(mae2)
#' sampleMap(mae2)
#'
#' @export renameColname
renameColname <- function(x, i, value) {
    stopifnot(length(i) == 1L, !is.na(i), !missing(i))
    exps <- experiments(x)
    expAssay <- x[[i]]
    splitmp <- mapToList(sampleMap(x))
    expmap <- splitmp[[i]]
    old <- expmap[["colname"]]
    if (length(old) != length(value))
        stop("'value' length does not match 'length(colnames(x[[i]]))'")
    newcns <- value[match(expmap[["colname"]], old)]
    expmap[["colname"]] <- newcns
    splitmp[[i]] <- expmap
    smp <- listToMap(splitmp)
    colnames(expAssay) <- value
    exps[[i]] <- expAssay
    BiocBaseUtils::setSlots(
        object = x,
        sampleMap = smp,
        ExperimentList = exps
    )
}

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases splitAssays
#'
#' @section splitAssays:
#' The `splitAssays` method separates columns in each of the assays based
#' on the `hitList` input. The `hitList` can be generated using
#' the `makeHitList` helper function. To use the `makeHitList`
#' helper, the user should input a list of patterns that will match on the
#' column names of each assay. These matches should be mutually exclusive as
#' to avoid repetition of columns across assays. See the examples section.
#'
#' @param hitList a named `list` or `List` of logical vectors that
#' indicate groupings in the assays
#'
#' @param patternList a named `list` or `List` of atomic character
#' vectors that are the input to `grepl` for identifying groupings in
#' the assays
#'
#' @export hasRowRanges
setGeneric("splitAssays", function(x, hitList)
        standardGeneric("splitAssays")
)

#' @rdname MultiAssayExperiment-helpers
#' @exportMethod splitAssays
setMethod("splitAssays", "MultiAssayExperiment",
    function(x, hitList) {
        stopifnot(is.list(hitList) || is(hitList, "List"))
        exps <- experiments(x)
        innames <- lapply(hitList, names)
        validENames <- unique(unlist(innames))
        exps <- exps[names(exps) %in% validENames]
        sublist <- lapply(hitList, function(logilist) {
            Map(function(x, y) {
                x[, y, drop = FALSE]
            }, x = exps, y = logilist)
        })
        sublist <- unlist(sublist, recursive = FALSE)
        names(sublist) <- gsub("\\.", "_", names(sublist))
        names(x) <- names(sublist)

        experiments(x) <- ExperimentList(sublist)
        x
    }
)

#' @rdname MultiAssayExperiment-helpers
#'
#' @aliases makeMatchList
#'
#' @examples
#'
#' patts <- list(
#'     normals = "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-11",
#'     tumors = "TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01"
#' )
#'
#' data("miniACC")
#'
#' hits <- makeHitList(miniACC, patts)
#'
#' ## only turmors present
#' splitAssays(miniACC, hits)
#'
#' @export
makeHitList <- function(x, patternList) {
    Colnames <- colnames(x)
    res <- lapply(
        .setNames(nm = names(patternList)),
        function(pattname, patt) {
            grepl(patt[[pattname]], unlist(Colnames))
        }, patt = patternList
    )
    reslist <- lapply(res, relist, Colnames)
    sums <- do.call(function(...) Map(`+`, ...), reslist)
    if (any(unlist(sums) > 1))
        stop("Groupings are not mutually exclusive")
    matching <- vapply(reslist, function(x) any(any(x)), logical(1L))
    reslist[matching]
}

