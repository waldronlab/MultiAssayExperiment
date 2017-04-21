#' @include RangedRaggedAssay-class.R MultiAssayExperiment-class.R
#' ExperimentList-class.R
#'
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom utils .DollarNames
#' @importFrom reshape2 melt
#' @importFrom tidyr gather
#' @importFrom IRanges SplitDataFrameList
NULL

.generateMap <- function(colData, experiments) {
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    matches <- match(colname, rownames(colData))
    if (length(matches) && all(is.na(matches)))
        stop("no way to map colData to ExperimentList")
    primary <- rownames(colData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)

    if (nrow(autoMap) && any(is.na(autoMap[["primary"]]))) {
        notFound <- autoMap[is.na(autoMap[["primary"]]), ]
        warning("Data from rows:",
                sprintf("\n %s - %s", notFound[, 2], notFound[, 3]),
                "\ndropped due to missing phenotype data")
        autoMap <- autoMap[!is.na(autoMap[["primary"]]), ]
    }
    autoMap
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get the dimension names for
#' an \code{ExperimentList} using \code{\link[IRanges]{CharacterList}}
setMethod("dimnames", "ExperimentList", function(x) {
    list(IRanges::CharacterList(lapply(x, rownames)),
    IRanges::CharacterList(lapply(x, colnames)))
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.checkFindOverlaps <- function(obj_cl) {
    return(
        all(hasMethod("findOverlaps", signature(obj_cl, "GRanges"),
                      where = c("package:GenomicRanges", "package:IRanges",
                                "package:SummarizedExperiment")),
            hasMethod("subsetByOverlaps", signature(obj_cl, "GRanges"),
                      where = c("package:GenomicRanges", "package:IRanges",
                                "package:SummarizedExperiment")))
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

#' @describeIn MultiAssayExperiment Subset a \code{MultiAssayExperiment} object
#' @param x A \code{MultiAssayExperiment} object for subsetting
#' @param i subsetting: Either a \code{character}, or \code{GRanges} object for subsetting
#' by rows, assay: unused argument (missing)
#' @param j Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by columns
#' @param k Either a \code{character}, \code{logical}, or \code{numeric} vector
#' for subsetting by assays
#' @param ... Additional arguments. See details for more information.
#' @param drop logical (default TRUE) whether to drop empty assay elements
#' in the \code{ExperimentList}
#' @aliases [,MultiAssayExperiment,ANY-method
setMethod("[", c("MultiAssayExperiment", "ANY", "ANY", "ANY"),
          .subsetMultiAssayExperiment)

#' @describeIn MultiAssayExperiment A \code{logical} value indicating an empty
#' \code{MultiAssayExperiment}
#' @exportMethod isEmpty
setMethod("isEmpty", "MultiAssayExperiment", function(x)
    length(x) == 0L)

#' Subset \code{MultiAssayExperiment} object by Assay type
#'
#' Select which assay(s) to obtain from available datasets
#'
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what assay(s) to select
#' @return A \code{\link{MultiAssayExperiment}} object
#' @seealso `subset,MultiAssayExperiment-method`
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
#' @export subsetByAssay
setGeneric("subsetByAssay", function(x, y) standardGeneric("subsetByAssay"))

#' @describeIn subsetByAssay Use either a \code{numeric}, \code{logical},
#' or \code{character} vector to subset assays in a
#' \code{MultiAssayExperiment}
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

#' @export
#' @rdname MultiAssayExperiment-methods
setMethod("[[", "MultiAssayExperiment", function(x, i, j, ...) {
    experiments(x)[[i]]
})

#' @export
#' @rdname MultiAssayExperiment-methods
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

#' Subset \code{MultiAssayExperiment} object by \code{colData} rows
#'
#' Select biological units in a \code{MultiAssayExperiment} with
#' \code{subsetByColData}
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what \code{colData} rows to select
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
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
#' @export subsetByColData
setGeneric("subsetByColData", function(x, y) standardGeneric("subsetByColData"))

#' @describeIn subsetByColData Either a \code{numeric}, \code{character}, or
#' \code{logical} vector to apply a column subset of a
#' \code{MultiAssayExperiment} object
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

#' @describeIn subsetByColData Use a \code{character}
#' vector for subsetting column names
setMethod("subsetByColData", c("MultiAssayExperiment", "character"),
          function(x, y) {
              y <- unique(y)
              if (!any(rownames(colData(x)) %in% y))
                  stop("No matching identifiers found")
              if (!all(y %in% rownames(colData(x))))
                  warning("Not all identifiers found in data")
              callNextMethod(x = x, y = y)
          })

#' Subset \code{MultiAssayExperiment} object
#'
#' \code{subsetByColumn} returns a subsetted
#' \code{\linkS4class{MultiAssayExperiment}} object
#'
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y Either a \code{numeric}, \code{character} or
#' \code{logical} object indicating what rownames in the colData to select
#' for subsetting
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' subsetByColumn(myMultiAssayExperiment, list(Affy = 1:2,
#'                 Methyl450k = c(3,5,2), RNASeqGene = 2:4, CNVgistic = 1))
#'
#' subsetWith <- mendoapply(`[`, colnames(myMultiAssayExperiment),
#'                         MoreArgs = list(1:2))
#' subsetByColumn(myMultiAssayExperiment, subsetWith)
#'
#' @export subsetByColumn
setGeneric("subsetByColumn", function(x, y) standardGeneric("subsetByColumn"))

#' @describeIn subsetByColumn Use a \code{list} to subset by
#' colname in a \code{MultiAssayExperiment}
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
    colData(x) <- colData(x)[rownames(colData(x)) %in% selectors, ]
    sampleMap(x) <- newMap
    return(x)
})

#' @describeIn subsetByColumn Use an S4 \code{List} to subset
#' a \code{MultiAssayExperiment}. The order of the subsetting
#' elements in this \code{List} must match that of the
#' \code{ExperimentList} in the \code{MultiAssayExperiment}.
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
            if (.checkFindOverlaps(class(element)))
                i <- IRanges::overlapsAny(element, i, ...)
            else
                i <- which(rownames(element) %in% as.character(i))
        } else if (is.character(i)) {
            i <- which(rownames(element) %in% i)
        } else if (!is.logical(i)) {
            i <- as.integer(i)
        }
        i
    }))
}

#' Subset \code{MultiAssayExperiment} object by Feature
#'
#' Subset a \code{MultiAssayExperiment} class by provided feature names or a
#' \code{GRanges} object
#'
#' @param x A \code{\link{MultiAssayExperiment}} object
#' @param y A \code{character} vector or \code{GRanges} class object
#' containing feature names or ranges
#' @param ... Additional arguments to pass to low level subsetting function
#' primarily when using a \code{GRanges} object for subsetting
#' (via \code{getHits})
#' @return A \code{\link{MultiAssayExperiment}} object
#'
#' @examples
#' ## Load a MultiAssayExperiment example
#' example("MultiAssayExperiment")
#'
#' ## Use a GRanges object to subset rows where ranged data present
#' egr <- GRanges(seqnames = "chr1", IRanges(start = 1, end = 3), strand = "-")
#' subsetByRow(myMultiAssayExperiment, egr)
#'
#' ## Use a logical vector (recycling used)
#' subsetByRow(myMultiAssayExperiment, c(TRUE, FALSE))
#'
#' ## Use a character vector
#' subsetByRow(myMultiAssayExperiment, "ENST00000355076")
#'
#' @export subsetByRow
setGeneric("subsetByRow", function(x, y, ...)
    standardGeneric("subsetByRow"))

#' @describeIn subsetByRow Subset a
#' \code{MultiAssayExperiment} with either a
#' \code{numeric} or \code{logical} vector
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

#' @describeIn subsetByRow Use a list of equal length as
#' the \code{ExperimentList} to subset. The order of
#' the subsetting elements in this list must match that of the
#' \code{ExperimentList} in the \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "list"), function(x, y) {
    if (length(x) != length(y))
        stop("List length must be the same as ExperimentList length")
    if (!identical(names(x), names(y)))
        stop("List input order much match that of the 'ExperimentList'")
    y <- as(y, "List")
    subsetByRow(x, y)
})

#' @describeIn subsetByRow Use an S4 \code{List} to subset
#' a \code{MultiAssayExperiment}. The order of the subsetting elements in
#' this \code{List} must match that of the \code{ExperimentList} in the
#' \code{MultiAssayExperiment}.
setMethod("subsetByRow", c("MultiAssayExperiment", "List"), function(x, y) {
    if (is(y, "DataFrame"))
        stop("Provide a list of indices for subsetting")
    if (is(y, "CharacterList"))
        y <- IRanges::IntegerList(mapply(function(expList, char) {
            which(rownames(expList) %in% char)
        }, expList = experiments(x), char = y))
    newExpList <- mendoapply(function(explist, i) {
        explist[i, , drop = FALSE]
    }, experiments(x), y)
    experiments(x) <- newExpList
    return(x)
})

#' @exportMethod complete.cases
#' @describeIn MultiAssayExperiment Return a logical vector of biological units
#' with data across all experiments
setMethod("complete.cases", "MultiAssayExperiment", function(...) {
    args <- list(...)
    if (length(args) == 1L) {
        oldMap <- sampleMap(args[[1L]])
        listMap <- mapToList(oldMap)
        allPrimary <- Reduce(intersect,
                             lapply(listMap,
                                    function(element) {
                                        element[["primary"]]
                                    }))
        rownames(colData(args[[1L]])) %in% allPrimary
    }
})

.splitArgs <- function(args) {
    assayArgNames <- c("mcolname", "background", "type",
                       "make.names", "ranges")
    assayArgs <- args[assayArgNames]
    altArgs <- args[!names(args) %in% assayArgNames]
    assayArgs <- Filter(function(x) !is.null(x), assayArgs)
    list(assayArgs, altArgs)
}

#' @describeIn MultiAssayExperiment Add an element to the
#' \code{ExperimentList} data slot
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
setMethod("c", "MultiAssayExperiment", function(x, ..., sampleMap = NULL,
                                                mapFrom = NULL) {
    newExperiments <- list(...)
    if (!length(newExperiments))
        stop("No arguments provided")
    if (is.list(newExperiments[[1L]]) || is(newExperiments[[1L]], "List") &&
        !is(newExperiments[[1L]], "DataFrame"))
        newExperiments <- ExperimentList(newExperiments[[1L]])
    else
        newExperiments <- ExperimentList(newExperiments)
    if (is.null(names(newExperiments)))
        stop("Additional experiments must be named")
    if (is.null(sampleMap)) {
        if (!is.null(mapFrom)) {
            warning("Assuming column order in the data provided ",
                    "\n matches the order in 'mapFrom' experiment(s) colnames")
            addMaps <- mapToList(sampleMap(x))[mapFrom]
            names(addMaps) <- names(newExperiments)
            sampleMap <- mapply(function(x, y) {
                x[["colname"]] <- colnames(y)
                return(x)
            }, addMaps, newExperiments)
        } else {
        sampleMap <- .generateMap(colData(x), newExperiments)
        }
    }
    if (is(sampleMap, "DataFrame") || is.data.frame(sampleMap))
        sampleMap <- mapToList(sampleMap)
    else if (!is.list(sampleMap))
        stop("Provided 'sampleMap' must be either a 'DataFrame' or a 'list'")
    newListMap <- c(mapToList(sampleMap(x)),
                    IRanges::SplitDataFrameList(sampleMap))
    sampleMap(x) <- listToMap(newListMap)
    experiments(x) <- c(experiments(x), newExperiments)
    validObject(x)
    return(x)
})
