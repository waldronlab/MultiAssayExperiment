#' @include MultiAssayExperiment-class.R ExperimentList-class.R
#'
#' @importFrom utils .DollarNames
NULL

.checkOverlapsAny <- function(obj_cl) {
    if (identical(obj_cl, c("matrix", "array")))
        obj_cl <- "matrix"
    return(any(
        hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
            getNamespace("GenomicRanges")),
        hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
            getNamespace("SummarizedExperiment")),
        hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
            getNamespace("IRanges"))
        )
    )
}

.rowIdx <- function(x) {
    IntegerList(lapply(x, function(exper) seq_len(dim(exper)[[1L]])))
}

.getHits <- function(expList, i, ...) {
    IntegerList(lapply(expList, function(element) {
        rnames <- rownames(element)
        if (is(i, "Vector")) {
            if (is(element, "RangedSummarizedExperiment"))
                element <- rowRanges(element)
            if (is(element, "VcfStack"))
                i <- which(rnames %in% as.character(seqnames(i)))
            if (.checkOverlapsAny(class(element)) &&
                !is(element, "SummarizedExperiment"))
                i <- which(overlapsAny(element, i, ...))
            else
                i <- match(intersect(as.character(i), rnames), rnames)
        } else if (is.character(i)) {
            i <- match(intersect(i, rnames), rnames)
        } else {
            i <- as.integer(i)
        }
        i
    }))
}

.matchReorderSub <- function(assayMap, identifiers) {
    positions <- unlist(
        lapply(
            identifiers,
            function(ident) {
                which(!is.na(match(assayMap[["primary"]], ident)))
            }
        )
    )
    assayMap[positions, ]
}

#' @name subsetBy
#'
#' @title Subsetting a MultiAssayExperiment object
#'
#' @description A set of functions for extracting and dividing a
#' `MultiAssayExperiment`
#'
#' @param x A `MultiAssayExperiment` or `ExperimentList`
#'
#' @param y Either a `character`, `integer`, `logical`, `list`, `List`,
#' or `GRanges` object for subsetting by rows _within the experiments_
#'
#' @param i For the `subsetByRow` and `subsetByRowData` `MultiAssayExperiment`
#' methods, either a `character`, `logical`, or `numeric` vector to selectively
#' subset experiments with `y` (default is `TRUE`). For **bracket** (`[`)
#' methods, see `y` input.
#'
#' @param j Either a `character`, `logical`, or `numeric` vector
#' for subsetting by `colData` rows. See details for more information.
#'
#' @param k Either a `character`, `logical`, or `numeric` vector
#' for subsetting by assays
#'
#' @param ... Additional arguments passed on to lower level functions.
#'
#' @param drop logical (default FALSE) whether to drop all empty assay elements
#' in the `ExperimentList`
#'
#' @param rowDataCol `character(1)` The name of the column in the `rowData`.
#' If the column is not present, the experiment will be skipped. When
#' `rowDataCol` is `"rownames"` or `"row.names"`, the values of `y` will
#' be matched with the row names in the `rowData` of the experiment.
#'
#' @aliases [,MultiAssayExperiment,ANY-method subsetByColData subsetByRow
#' subsetByColumn subsetByAssay subset subsetBy
#'
#' @details
#' Subsetting a MultiAssayExperiment by the **j** index can yield a call
#' to either `subsetByColData` or `subsetByColumn`. For vector inputs,
#' the subset will be applied to the `colData` rows. For `List`-type
#' inputs, the List will be applied to each of the elements in the
#' `ExperimentList`.
#' The order of the subsetting elements in the
#' `List` must match that of the `ExperimentList` in the
#' `MultiAssayExperiment`.
#'
#' * `subsetBycolData`: Select biological units by vector input types
#' * `subsetByColumn`: Select observations by assay or for each assay
#' * `subsetByRow`: Select rows by assay or for each assay
#' * `subsetByAssay`: Select experiments
#'
#' @return `subsetBy*`: operations are endomorphic and return either
#' `MultiAssayExperiment` or `ExperimentList` depending on the
#' input.
#'
#' @examples
#' ## Load the example MultiAssayExperiment
#' example("MultiAssayExperiment")
#'
#' ## Using experiment names
#' subsetByAssay(mae, "Affy")
#'
#' ## Using numeric indices
#' subsetByAssay(mae, 1:2)
#'
#' ## Using a logical vector
#' subsetByAssay(mae, c(TRUE, FALSE, TRUE))
#'
#' ## Subset by character vector (Jack)
#' subsetByColData(mae, "Jack")
#'
#' ## Subset by numeric index of colData rows (Jack and Bob)
#' subsetByColData(mae, c(1, 3))
#'
#' ## Subset by logical indicator of colData rows (Jack and Jill)
#' subsetByColData(mae, c(TRUE, TRUE, FALSE, FALSE))
#'
#' subsetByColumn(mae, list(Affy = 1:2,
#'     Methyl450k = c(3,5,2), RNASeqGene = 2:4, GISTIC = 1))
#'
#' subsetWith <- S4Vectors::mendoapply(`[`, colnames(mae),
#'     MoreArgs = list(1:2))
#' subsetByColumn(mae, subsetWith)
#'
#' ## Use a GRanges object to subset rows where ranged data present
#' egr <- GenomicRanges::GRanges(seqnames = "chr2",
#'     IRanges::IRanges(start = 11, end = 13), strand = "-")
#' subsetByRow(mae, egr)
#'
#' ## Use a logical vector (recycling used)
#' subsetByRow(mae, c(TRUE, FALSE))
#'
#' ## Use a character vector
#' subsetByRow(mae, "ENST00000355076")
#'
#' ## Use i index to selectively subsetByRow
#' subsetByRow(mae, "ENST00000355076", i = c(TRUE, TRUE, FALSE, FALSE))
#'
#' ## Use i index to selectively subsetByRowData
#' subsetByRowData(
#'     mae, "ENST00000355076", "rownames", i = "Affy"
#' )
NULL

# subsetBy Generics -------------------------------------------------------

#' @rdname subsetBy
#' @export subsetByRow
setGeneric("subsetByRow", function(x, y, ...) standardGeneric("subsetByRow"))

#' @rdname subsetBy
#' @export subsetByRowData
setGeneric(
    "subsetByRowData",
    function(x, y, rowDataCol, ...) standardGeneric("subsetByRowData")
)

#' @rdname subsetBy
#' @export subsetByColData
setGeneric("subsetByColData", function(x, y) standardGeneric("subsetByColData"))

#' @rdname subsetBy
#' @export subsetByColumn
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))

#' @rdname subsetBy
#' @export subsetByAssay
setGeneric("subsetByAssay", function(x, y) standardGeneric("subsetByAssay"))

.subsetCOLS <- function(object, cutter) {
    mendoapply(function(x, j) {
        if (!is.null(j))
            x[, j, drop = FALSE]
        else
            x
    }, x = object, j = cutter)
}

.subsetROWS <- function(object, cutter) {
    mendoapply(function(x, i) {
        if (!is.null(rownames(x)) && !is.null(i))
            x[i, , drop = FALSE]
        else
            x
    }, x = object, i = cutter)
}

.fillEmptyExps <- function(exps, subr) {
    if (!any(names(subr) %in% names(exps)))
        stop("No matching experiment names in subset list", call. = FALSE)
    if (!all(names(exps) %in% names(subr))) {
        outnames <- setdiff(names(exps), names(subr))
        names(outnames) <- outnames
        subr <- c(subr, lapply(outnames, function(x) NULL))
    }
    subr[names(exps)]
}

# subsetByRow,ExperimentList-methods -----------------------------------------

#' @rdname subsetBy
setMethod("subsetByRow", c("ExperimentList", "ANY"), function(x, y, ...) {
    rowIds <- .rowIdx(x)
    if (is.integer(y)) {
        lowerLimit <- min(max(rowIds))
        if (max(y) > lowerLimit)
            stop("subscript contains out-of-bounds indices,\n",
                " use an ", sQuote("IntegerList"), " index for finer control")
    }
    subsetor <- .getHits(x, y, ...)
    subsetByRow(x, subsetor)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("ExperimentList", "list"), function(x, y, ...) {
    y <- .fillEmptyExps(x, y)
    .subsetROWS(x, y)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("ExperimentList", "List"), function(x, y, ...) {
    if (is(y, "DataFrame") || is(y, "GRangesList"))
        stop("Provide a list of indices for subsetting")
    if (is(y, "GRanges"))
        return(callNextMethod())
    y <- as.list(y)
    subsetByRow(x, y)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("ExperimentList", "logical"), function(x, y, ...) {
    logi <- LogicalList(rep(list(y), length(x)))
    x[logi]
})

# subsetByColumn,ExperimentList-methods -----------------------------------

#' @rdname subsetBy
setMethod("subsetByColumn", c("ExperimentList", "list"), function(x, y) {
    y <- .fillEmptyExps(x, y)
    .subsetCOLS(x, y)
})

#' @rdname subsetBy
setMethod("subsetByColumn", c("ExperimentList", "List"), function(x, y) {
    Y <- as.list(y)
    subsetByColumn(x, Y)
})

#' @rdname subsetBy
setMethod("subsetByColumn", c("ExperimentList", "logical"), function(x, y) {
    Y <- endoapply(colnames(x), `[`, y)
    .subsetCOLS(x, Y)
})

# subsetByAssay,ExperimentList-methods ------------------------------------

#' @rdname subsetBy
setMethod("subsetByAssay", c("ExperimentList", "ANY"), function(x, y) {
    x[y]
})

# subsetByColData,MultiAssayExperiment-methods ----------------------------

#' @rdname subsetBy
setMethod("subsetByColData", c("MultiAssayExperiment", "ANY"), function(x, y) {
    coldata <- colData(x)
    if (length(y) > nrow(coldata))
        stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
            call. = FALSE)
    newcoldata <- coldata[y, , drop = FALSE]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(
        listMap,
        function(elementMap, keepers) {
            .matchReorderSub(
                elementMap, intersect(keepers, elementMap[["primary"]])
            )
        }, keepers = rownames(newcoldata)
    )
    newMap <- listToMap(listMap, fill = FALSE)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    columns <- .fillEmptyExps(experiments(x), columns)
    newSubset <- Map(
        function(x, j) {
            x[, j, drop = FALSE]
        },
        x = experiments(x), j = columns
    )
    newSubset <- ExperimentList(newSubset)

    BiocBaseUtils::setSlots(x,
        ExperimentList = newSubset,
        colData = newcoldata,
        sampleMap = newMap,
        check = FALSE
    )
})

#' @rdname subsetBy
setMethod("subsetByColData", c("MultiAssayExperiment", "character"),
    function(x, y) {
        coldata <- colData(x)
        if (length(y) > nrow(coldata))
            stop("subscript vector 'j' in 'mae[i, j, k]' is out-of-bounds",
                call. = FALSE)
        y <- unique(y)
        if (!any(rownames(colData(x)) %in% y))
            stop("No matching 'colData' row identifiers provided")
        if (!all(y %in% rownames(colData(x))))
            warning("Not all provided identifiers found in 'colData'")
        callNextMethod(x = x, y = y)
})

# subsetByRow,MultiAssayExperiment-method ---------------------------------

#' @rdname subsetBy
#' @exportMethod subsetByRow
setMethod(
    "subsetByRow", c("MultiAssayExperiment", "ANY"),
    function(x, y, i = TRUE, ...) {
        stopifnot(
            !anyNA(i) &&
                (is.logical(i) || is.character(i) || is.numeric(i))
        )
        experiments(x)[i] <- subsetByRow(experiments(x)[i], y)
        return(x)
    }
)

#' @rdname subsetBy
setMethod(
    "subsetByRow", c("MultiAssayExperiment", "list"),
    function(x, y, ...) {
        experiments(x) <- subsetByRow(experiments(x), y)
        return(x)
    }
)

#' @rdname subsetBy
setMethod(
    "subsetByRow", c("MultiAssayExperiment", "List"),
    function(x, y, ...) {
        y <- as.list(y)
        subsetByRow(x, y)
    }
)

# subsetByColumn,MultiAssayExperiment-method ------------------------------

#' @rdname subsetBy
setMethod("subsetByColumn", c("MultiAssayExperiment", "ANY"), function(x, y) {
    if (is.character(y) || is.logical(y) || is.numeric(y))
        subsetByColData(x, y)
    else {
        experiments(x) <- subsetByColumn(experiments(x), y)
        return(x)
    }
})

# subsetByAssay,MultiAssayExperiment-method -------------------------------

#' @rdname subsetBy
setMethod("subsetByAssay", c("MultiAssayExperiment", "ANY"), function(x, y) {
    subexp <- experiments(x)[y]
    dropnames <- setdiff(names(experiments(x)), names(subexp))
    if (length(dropnames)) {
        if (isEmpty(drops(x)))
            warning("'experiments' dropped; see 'drops()'", call. = FALSE)
        drops(x) <- list(experiments = dropnames)
    }
    experiments(x) <- subexp
    return(x)
})

# subsetByRowData,MultiAssayExperiment-method -----------------------------

#' @rdname subsetBy
setMethod(
    "subsetByRowData", c("MultiAssayExperiment", "character", "character"),
    function(x, y, rowDataCol, i = TRUE, ...) {
        if (is.character(i))
            logi <- names(x) %in% i
        else if (is.logical(i) || is.numeric(i))
            logi <- names(x) %in% names(x)[i]
        else
            stop("Invalid experiment subscript type for 'i'")
        valids <- hasRowData(x)[which(logi)]
        if (any(!valids)) {
            notValids <- paste(
                names(valids[!valids]), collapse = ", "
            )
            stop("Selected experiments have no 'rowData': ", notValids)
        }
        i <- hasRowData(x) & logi
        if (!any(i))
            stop("No 'rowData' available for subsetting")
        y <- lapply(
            experiments(x)[i],
            function(exper) {
                rd <- rowData(exper)
                if (rowDataCol %in% c("rownames", "row.names"))
                    rownames(rd) %in% y
                else if (rowDataCol %in% colnames(rd))
                    rd[[rowDataCol]] %in% y
                else
                    TRUE
            }
        )
        subsetByRow(x = x, y = y, i = i)
    }
)
