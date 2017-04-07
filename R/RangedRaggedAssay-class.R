### ==============================================
### RangedRaggedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#'
#' @exportClass RangedRaggedAssay
#' @name RangedRaggedAssay-class
#'
#' @docType class
NULL

#' @keywords internal
.RangedRaggedAssay <- setClass("RangedRaggedAssay", contains = "GRangesList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Create a RangedRaggedAssay
#'
#' Construct an object representing ranged-based data, typically from a
#' \code{\link{GRangesList}}. The assay method will extract a particular column
#' from the metadata and represent it in a matrix. See the \code{show} method
#' for an example.
#'
#' @section Deprecated:
#' The \code{RangedRaggedAssay} class is \strong{deprecated} and defunct by the next
#' release cycle. Please use the \strong{RaggedExperiment} class to represent
#' copy number, mutation and other genomic range based data. See
#' \code{RaggedExperiment} for more detail.
#'
#' @param x A \code{list}, \code{GRanges} or \code{GRangesList} object
#' @return A \code{\linkS4class{RangedRaggedAssay}} class object
#'
#' @example inst/scripts/RangedRaggedAssay-Ex.R
#'
#' @seealso \code{\link{assay,RangedRaggedAssay,missing-method}}
#'
#' @export RangedRaggedAssay
RangedRaggedAssay <- function(x = GRangesList()) {
    .Deprecated("RaggedExperiment")
    if (is(x, "GRanges")) {
        x <- GRangesList(x)
    }
    if (is(x, "GRangesList")) {
        metad <- mcols(x)
        missingRownames <- vapply(X = x, FUN = function(grl) {
            is.null(names(grl))
        }, FUN.VALUE = logical(1L))
        if (all(missingRownames)) {
            u_obj <- unlist(x, use.names = FALSE)
            gre <- granges(u_obj, use.mcols = TRUE)
            names(gre) <- as.character(gre)
            x <- relist(gre, x)
        }
    }
    newRRA <- .RangedRaggedAssay(x)
    mcols(newRRA) <- metad
    return(newRRA)
}


.RangedBracketSubsetRRA <- function(x, i, j, ..., drop) {
    if (length(drop) != 1L || (!missing(drop) && drop)) {
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
    }
    if (!missing(j)) {
        x <- callNextMethod(x = x, i = j)
    }
    if (!missing(i)) {
        x <- endoapply(x, function(rra) {
            IRanges::subsetByOverlaps(rra, i, ...)
            # x <- x[relist(subsetByOverlaps(unlist(x,
            # use.names = FALSE), i, ...), x)]
        })
    }
    return(x)
}

.sBracketSubsetRRA <- function(x, i, j, ..., drop) {
    if (length(drop) != 1L || (!missing(drop) && drop)) {
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")
    }
    if (missing(i) && missing(j)) {
        return(x)
    }
    if (!missing(j)) {
        x <- callNextMethod(x = x, i = j)
    }
    if (!missing(i)) {
        if (is.character(i)) {
            cLL <- relist(names(unlist(x, use.names = FALSE)) %in% i, x)
            x <- callNextMethod(x = x, i = cLL)
        } else if (is.numeric(i) || is.logical(i)) {
            sampleFactor <- rep(names(x), lengths(x))[i]
            x <- unlist(x, use.names = FALSE)[i]
            x <- RangedRaggedAssay(IRanges::splitAsList(x, sampleFactor))
        } else {
            x <- callNextMethod(x = x, i = i)
        }
    }
    return(x)
}

#' Subset RangedRaggedAssay
#'
#' @description
#' Subsetting a RangedRaggedAssay can be done using either rownames and column
#' names
#'
#' @details
#' The \code{\ldots} argument allows the user to specify arguments in the
#' \code{\link{subsetByOverlaps}} function. When calling the \code{reduce}
#' method, the additional arguments correspond to those in either the
#' \code{assay} method or the \code{reduce} method. The \code{reduce} arguments
#' include a function for applying over the rows (\code{combine}) and a
#' \code{vectorized} argument which indicates whether the given function is
#' vectorized or not.
#'
#' @param x A \code{\link{RangedRaggedAssay}} class
#' @param i Either a \code{character} or \code{GRanges} class object
#' to subset by rows
#' @param j Either a \code{character}, \code{numeric}, or \code{logical}
#' type for selecting columns (\code{\link[GenomicRanges]{GRangesList}} method)
#' @param ... Additional arguments. See details for more information.
#' @param drop logical (default TRUE) whether to drop empty columns
#' @seealso \code{\link[IRanges]{findOverlaps-methods}}
#' @return A \code{\link{RangedRaggedAssay}} class object
#' @describeIn RangedRaggedAssay Subset a \code{RangedRaggedAssay} with either
#' \code{chracter}, \code{numeric}, or \code{logical}
#' @aliases [,RangedRaggedAssay,ANY-method
setMethod("[", c("RangedRaggedAssay", "ANY", "ANY"),
          .sBracketSubsetRRA)

#' @describeIn RangedRaggedAssay Subset a \code{RangedRaggedAssay} using a
#' \code{GRanges} class object
#' @aliases [,RangedRaggedAssay,GRanges-method
setMethod("[", c("RangedRaggedAssay", "GRanges", "ANY"),
          .RangedBracketSubsetRRA)

#' @describeIn RangedRaggedAssay Obtain dimension lengths of a
#' \code{RangedRaggedAssay} class object
setMethod("dim", "RangedRaggedAssay", function(x)
    c(length(unlist(x)), length(x)))

#' @describeIn RangedRaggedAssay Get dimension names
#' for a \code{RangedRaggedAssay}
setMethod("dimnames", "RangedRaggedAssay", function(x) {
    dgr <- names(unlist(x, use.names = FALSE))
    list(dgr, names(x))
})

#' @exportMethod dimnames<-
#' @describeIn RangedRaggedAssay value: A modified \code{RangedRaggedAssay}
#' object
#' @param value A \code{list} object of row and column names
setReplaceMethod("dimnames", c("RangedRaggedAssay", "list"),
    function(x, value)
{
    names(x@unlistData) <- value[[1]]
    names(x) <- value[[2]]
    x
})

.findNumericMcol <- function(rra) {
    mcolnames <- names(mcols(rra[[1L]]))
    meanCols <- grepl("mean", mcolnames, ignore.case = TRUE)
    scoreCol <- grepl("score", mcolnames, ignore.case = TRUE)
    if (any(meanCols)) {
        if (length(mcolnames[meanCols]) > 1L)
            warning("More than one 'mcol' found,",
                    " selecting the first")
        nummcol <- mcolnames[meanCols][[1L]]
    } else if (any(scoreCol)) {
        nummcol <- mcolnames[scoreCol]
        stopifnot(S4Vectors::isSingleString(nummcol))
    } else {
        stop("No numeric 'mcolname' specified")
    }
    return(nummcol)
}

#' @describeIn RangedRaggedAssay Separate non-disjoint ranges and merge
#' with function
#' @param mcolname A single character string indicating metadata column to use
#' for summaries
#' @importFrom IRanges disjoin
#' @exportMethod disjoin
setMethod("disjoin", "RangedRaggedAssay",
          function(x, mcolname = NULL, simplify = BiocGenerics::mean, ...) {
    .Deprecated("RaggedExperiment")
    if (is.null(mcolname))
        mcolname <- .findNumericMcol(x)
    if (any(!isDisjoint(x))) {
    newX <- lapply(x, function(singleRange, summarizer) {
        mCols <- mcols(singleRange)
        mCols <- mCols[, -which(names(mCols) == mcolname),
                       drop = FALSE]
        dj <- disjoin(singleRange, with.revmap=TRUE)
        revMap <- mcols(dj)
        revMap[, mcolname] <- sapply(dj$revmap, function(i) {
            summaryScore <- summarizer(mcols(singleRange)[[mcolname]][i])
            summaryScore
        })
        mcols(dj) <- revMap
        return(dj)
    }, summarizer = simplify)
    return(RangedRaggedAssay(GRangesList(newX)))
    }
    return(x)
})

#' @exportMethod show
#' @describeIn RangedRaggedAssay show method for
#' the \code{RangedRaggedAssay} class
#' @param object A \code{RangedRaggedAssay} class object
setMethod("show", "RangedRaggedAssay", function(object) {
    # if (!all(GenomicRanges::isDisjoint(object))) {
    #     cat("Non-disjoint RangedRaggedAssay")
    #     callNextMethod(object)
    # }
    metacols <- mcols(unlist(object))
    showable <- vapply(metacols, function(mcol) {
        is.atomic(mcol)
    }, logical(1L))
    elts <- names(metacols)
    elts <- elts[showable]

    cat(class(object), "with",
        length(dimnames(object)[[1]]), "disjoint ranges,",
        length(dimnames(object)[[2]]), "samples, and",
        length(elts), "data elements")
    # length(elts), "data elements\n")

    # for (elt in head(elts, 3)) {
    #     cat("\n", elt, "\n", sep="")
    #     x <- assay(object, mcolname = elt)
    #     print(head(x, 3))
    #     if (nrow(x) > 3)
    #         cat("...\n")
    # }
    # if (length(elts) > 3)
    #     cat("\n...")
    cat("\n")
})
