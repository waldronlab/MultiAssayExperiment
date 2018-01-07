#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' ExperimentList-class.R
#'
#' @import BiocGenerics S4Vectors methods
#' @importFrom BiocGenerics duplicated
#' @importFrom utils .DollarNames
#' @importFrom stats kmeans
#' @importFrom reshape2 melt
#' @importFrom tidyr gather
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges findOverlaps subsetByOverlaps
#' @importFrom IRanges splitAsList SplitDataFrameList endoapply mendoapply
#' @importFrom S4Vectors endoapply mendoapply
#' @importFrom IRanges IntegerList CharacterList LogicalList
#' @importFrom SummarizedExperiment findOverlaps assays
NULL

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.checkOverlapsAny <- function(obj_cl) {
    return(
        any(hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
                      getNamespace("GenomicRanges")),
            hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
                      getNamespace("SummarizedExperiment")),
            hasMethod("overlapsAny", signature(obj_cl, "GRanges"),
                      getNamespace("IRanges")))
    )
}

.isEmpty <- function(object) {
    isTRUE(unname(dim(object)[1]) == 0L || unname(dim(object)[2]) == 0L)
}

.subsetMultiAssayExperiment <- function(x, i, j, k, ..., drop = TRUE) {
    if (missing(i) && missing(j) && missing(k)) {
        return(x)
    }
    if (!missing(k)) {
        x <- subsetByAssay(x, k)
    }
    if (!missing(j)) {
        if (is(j, "list") || is(j, "List"))
            x <- subsetByColumn(x, j)
        else
            x <- subsetByColData(x, j)
    }
    if (!missing(i)) {
        x <- subsetByRow(x, i, ...)
    }
    if (drop) {
        isEmptyAssay <- vapply(experiments(x), FUN = .isEmpty,
                               FUN.VALUE = logical(1L))
        if (all(isEmptyAssay)) {
            experiments(x) <- ExperimentList()
        } else if (any(isEmptyAssay)) {
            keeps <- names(isEmptyAssay)[
                vapply(isEmptyAssay, function(k) {
                    !isTRUE(k)}, logical(1L))]
            x <- x[ , , keeps, drop = TRUE]
        }
    }
    return(x)
}

#' @name subsetBy
#' @title Subsetting a MultiAssayExperiment object
#' @description A set of functions for extracting and dividing a
#' \code{MultiAssayExperiment}
#'
#' @param x A \code{MultiAssayExperiment}
#' @param i Either a \code{character}, \code{integer}, \code{logical} or
#' \code{GRanges} object for subsetting by rows
#' @param j Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by \code{colData} rows. See details for more information.
#' @param k Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by assays
#' @param ... Additional arguments passed on to lower level functions.
#' @param drop logical (default TRUE) whether to drop empty assay elements
#' in the \code{ExperimentList}
#'
#' @aliases [,MultiAssayExperiment,ANY-method subsetByColData subsetByRows
#' subsetByColumns subsetByAssay subset subsetBy
#'
#' @details
#' Subsetting a MultiAssayExperiment by the \strong{j} index can yield a call
#' to either \code{subsetByColData} or \code{subsetByColumn}. For vector inputs,
#' the subset will be applied to the \code{colData} rows. For \code{List}-type
#' inputs, the List will be applied to each of the elements in the
#' \code{ExperimentList}.
#' The order of the subsetting elements in the
#' \code{List} must match that of the \code{ExperimentList} in the
#' \code{MultiAssayExperiment}.
#'
#' \itemize{
#' \item subsetBycolData: Select biological units by vector input types
#' \item subsetByColumn: Select observations by assay or for each assay
#' \item subsetByRow: Select rows by assay or for each assay
#' \item subsetByAssay: Select experiments
#' }
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Using experiment names
#' subsetByAssay(myMultiAssayExperiment, "Affy")
#'
#' ## Using numeric indicators
#' subsetByAssay(myMultiAssayExperiment, 1:2)
#'
#' ## Using a logical vector
#' subsetByAssay(myMultiAssayExperiment, c(TRUE, FALSE, TRUE))
#'
#' ## Subset by character vector (Jack)
#' subsetByColData(myMultiAssayExperiment, "Jack")
#'
#' ## Subset by numeric index of colData rows (Jack and Bob)
#' subsetByColData(myMultiAssayExperiment, c(1, 3))
#'
#' ## Subset by logical indicator of colData rows (Jack and Jill)
#' subsetByColData(myMultiAssayExperiment, c(TRUE, TRUE, FALSE, FALSE))
#'
#' subsetByColumn(myMultiAssayExperiment, list(Affy = 1:2,
#'                 Methyl450k = c(3,5,2), RNASeqGene = 2:4, CNVgistic = 1))
#'
#' subsetWith <- S4Vectors::mendoapply(`[`, colnames(myMultiAssayExperiment),
#'                         MoreArgs = list(1:2))
#' subsetByColumn(myMultiAssayExperiment, subsetWith)
#'
#' ## Use a GRanges object to subset rows where ranged data present
#' egr <- GenomicRanges::GRanges(seqnames = "chr1",
#'     IRanges::IRanges(start = 1, end = 3), strand = "-")
#' subsetByRow(myMultiAssayExperiment, egr)
#'
#' ## Use a logical vector (recycling used)
#' subsetByRow(myMultiAssayExperiment, c(TRUE, FALSE))
#'
#' ## Use a character vector
#' subsetByRow(myMultiAssayExperiment, "ENST00000355076")
#'
NULL

#' @rdname subsetBy
#' @aliases [,MultiAssayExperiment,ANY,ANY,ANY-method
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
          .subsetMultiAssayExperiment)

#' @export
#' @rdname subsetBy
setMethod("[[", "MultiAssayExperiment", function(x, i, j, ...) {
    experiments(x)[[i]]
})

#' @export
#' @rdname subsetBy
#' @param value An assay compatible with the MultiAssayExperiment API
setReplaceMethod("[[", "MultiAssayExperiment", function(x, i, j, ..., value) {
                         if (!missing(j) || length(list(...)) > 0)
                             stop("invalid replacement")
                         origLen <- length(x)
                         experiments(x) <- S4Vectors::setListElement(
                             experiments(x),
                             i, value)
                         if (origLen < length(x))
                            stop("replacement length greater than original")
                         return(x)
})

.matchReorderSub <- function(assayMap, identifiers) {
    positions <- unlist(
        lapply(identifiers,
               function(ident) {
                   which(!is.na(match(assayMap[["primary"]], ident)))
               }))
    assayMap[positions, ]
}

#' @rdname subsetBy
#' @export subsetByColData
setGeneric("subsetByColData", function(x, y) standardGeneric("subsetByColData"))

#' @rdname subsetBy
setMethod("subsetByColData", c("MultiAssayExperiment", "ANY"), function(x, y) {
    if (is.logical(y) || is.numeric(y))
        y <- unique(rownames(colData(x))[y])
    selectors <- y[y %in% rownames(colData(x))]
    newcolData <- colData(x)[
        match(selectors, rownames(colData(x))), , drop = FALSE]
    listMap <- mapToList(sampleMap(x), "assay")
    listMap <- lapply(listMap, function(elementMap, keepers) {
        .matchReorderSub(elementMap, keepers)
    }, keepers = selectors)
    newMap <- listToMap(listMap)
    columns <- lapply(listMap, function(mapChunk) {
        mapChunk[, "colname", drop = TRUE]
    })
    newSubset <- mapply(function(x, j) {x[, j, drop = FALSE]},
                        x = experiments(x), j = columns, SIMPLIFY = FALSE)
    newSubset <- ExperimentList(newSubset)
    experiments(x) <- newSubset
    sampleMap(x) <- newMap
    colData(x) <- newcolData
    return(x)
})

#' @rdname subsetBy
setMethod("subsetByColData", c("MultiAssayExperiment", "character"),
          function(x, y) {
              y <- unique(y)
              if (!any(rownames(colData(x)) %in% y))
                  stop("No matching identifiers found")
              if (!all(y %in% rownames(colData(x))))
                  warning("Not all identifiers found in data")
              callNextMethod(x = x, y = y)
          })

#' @rdname subsetBy
#' @export subsetByColumn
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))

#' @rdname subsetBy
setMethod("subsetByColumn", c("MultiAssayExperiment", "list"), function(x, y)
{
    y <- y[names(x)]
    experiments(x) <- ExperimentList(mapply(function(expList, j) {
        expList[, j, drop = FALSE]
    }, expList = experiments(x), j = y, SIMPLIFY = FALSE))
    newSamps <- as.list(colnames(x))
    listMap <- mapToList(sampleMap(x), "assay")
    newMap <- mapply(function(lMap, nSamps) {
        lMap[na.omit(match(nSamps,
                           as.character(lMap[["colname"]]))), ]
    }, lMap = listMap, nSamps = newSamps, SIMPLIFY = FALSE)
    newMap <- listToMap(newMap)
    selectors <- unique(as.character(newMap[["primary"]]))
    colData(x) <- colData(x)[rownames(colData(x)) %in% selectors, , drop = FALSE]
    sampleMap(x) <- newMap
    return(x)
})

#' @rdname subsetBy
setMethod("subsetByColumn", c("MultiAssayExperiment", "List"),
          function(x, y) {
              Y <- as.list(y)
              subsetByColumn(x, Y)
          })

.rowIdx <- function(x) {
    IRanges::IntegerList(lapply(x, function(exper) seq_len(dim(exper)[[1L]])))
}

.getHits <- function(expList, i, ...) {
    IRanges::IntegerList(lapply(expList, function(element) {
        if (is(i, "Vector")) {
            if (is(element, "RangedSummarizedExperiment"))
                element <- rowRanges(element)
            if (is(element, "VcfStack"))
                i <- which(rownames(element) %in% as.character(seqnames(i)))
            if (.checkOverlapsAny(class(element)) &&
                !is(element, "SummarizedExperiment"))
                i <- which(IRanges::overlapsAny(element, i, ...))
            else
                i <- match(intersect(as.character(i), rownames(element)),
                           rownames(element))
                # i <- na.omit(match(rownames(element), as.character(i)))
        } else if (is.character(i)) {
            i <- match(intersect(i, rownames(element)), rownames(element))
            # i <- na.omit(match(rownames(element), i))
        } else {
            i <- as.integer(i)
        }
        i
    }))
}

#' @rdname subsetBy
#' @export subsetByRow
setGeneric("subsetByRow", function(x, y, ...)
    standardGeneric("subsetByRow"))

#' @rdname subsetBy
setMethod("subsetByRow", c("MultiAssayExperiment", "ANY"), function(x, y, ...) {
    rowIds <- .rowIdx(experiments(x))
    if (is.integer(y)) {
        lowerLimit <- min(max(rowIds))
        if (max(y) > lowerLimit)
            stop("subscript contains out-of-bounds indices,\n",
                 " use an ", sQuote("IntegerList"), " index for finer control")
    }
    subsetor <- .getHits(experiments(x), y, ...)
    y <- rowIds[subsetor]
    subsetByRow(x, y)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("MultiAssayExperiment", "logical"), function(x, y) {
    rowIds <- .rowIdx(experiments(x))
    y <- S4Vectors::endoapply(rowIds, `[`, y)
    subsetByRow(x, y)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("MultiAssayExperiment", "list"), function(x, y) {
    if (length(x) != length(y))
        stop("List length must be the same as ExperimentList length")
    if (!identical(names(x), names(y)))
        stop("List input order much match that of the 'ExperimentList'")
    y <- as(y, "List")
    subsetByRow(x, y)
})

#' @rdname subsetBy
setMethod("subsetByRow", c("MultiAssayExperiment", "List"), function(x, y) {
    if (is(y, "DataFrame"))
        stop("Provide a list of indices for subsetting")
    if (is(y, "CharacterList"))
        y <- IRanges::IntegerList(Map(function(expList, char) {
            which(rownames(expList) %in% char)
        }, expList = experiments(x), char = y))
    newExpList <- S4Vectors::mendoapply(function(explist, i) {
        explist[i, , drop = FALSE]
    }, experiments(x), y)
    experiments(x) <- newExpList
    return(x)
})

#' @rdname subsetBy
#' @export subsetByAssay
#' @param y Any argument used for subsetting, can be a \code{character},
#' \code{logical}, \code{integer}, \code{list} or \code{List} vector
setGeneric("subsetByAssay", function(x, y) standardGeneric("subsetByAssay"))

#' @rdname subsetBy
setMethod("subsetByAssay", c("MultiAssayExperiment", "ANY"), function(x, y) {
    newSubset <- experiments(x)[y]
    listMap <- mapToList(sampleMap(x), "assay")
    ## TODO: Add sensible error message here
    newMap <- listMap[y]
    newMap <- listToMap(newMap)
    sampleMap(x) <- newMap
    experiments(x) <- newSubset
    return(x)
})

