#' @include MultiAssayExperiment-class.R
NULL

### ==============================================
### MatchedAssayExperiment class
### ----------------------------------------------

#' MatchedAssayExperiment - A matched-samples MultiAssayExperiment class
#'
#' @description
#' This class supports the use of matched samples where an equal number
#' of observations per biological unit are present in all assays.
#'
#' @return A \code{MatchedAssayExperiment} object
#'
#' @exportClass MatchedAssayExperiment
#' @seealso \link{MultiAssayExperiment}
#'
setClass("MatchedAssayExperiment", contains="MultiAssayExperiment")

.checkEqualPrimaries <- function(object) {
    listMap <- mapToList(sampleMap(object))
    primaryIDs <- lapply(listMap, function(x) x[["primary"]])
    allIDsEqual <- all(vapply(seq_along(primaryIDs)[-1], function(i, prim) {
        identical(prim[[1L]], prim[[i]])
    }, FUN.VALUE = logical(1L), prim = primaryIDs))
    if (!allIDsEqual)
        "Primary identifiers are not equal across assays"
    else
        NULL
}

.checkPrimaryOrder <- function(object) {
    colPrimary <- rownames(colData(object))
    listMap <- mapToList(sampleMap(object))
    primaryIDs <- lapply(listMap, function(x) x[["primary"]])
    allOrdered <- all(vapply(primaryIDs, function(prim) {
        identical(colPrimary, prim)
    }, logical(1L)))
    if (!allOrdered)
        "colData row identifiers not identical to sampleMap primary column"
    else
        NULL
}

.validMatchedAssayExperiment <- function(object) {
    if (length(object) != 0L) {
        c(.checkEqualPrimaries(object), .checkPrimaryOrder(object))
    }
}

S4Vectors::setValidity2("MatchedAssayExperiment", .validMatchedAssayExperiment)

.doMatching <- function(from) {
    if (!isEmpty(from)) {
        from <- intersectColumns(from)

        if (all(!lengths(colnames(from))))
            stop("No biological unit(s) measured across all assays")

        if (any(anyReplicated(from)))
            stop("Resolve replicate columns")
    }
    from
}

#' @describeIn MatchedAssayExperiment-class Construct a
#' \code{MatchedAssayExperiment} class from \linkS4class{MultiAssayExperiment}
#'
#' @param ... Either a single MultiAssayExperiment or the components to create
#' a valid MultiAssayExperiment
#'
#' @examples
#' data("miniACC")
#' acc <- as(miniACC, "MatchedAssayExperiment")
#' acc
#'
#' @aliases MatchedAssayExperiment
#' coerce,MultiAssayExperiment,MatchedAssayExperiment-method
#'
#' @export MatchedAssayExperiment
MatchedAssayExperiment <- function(...) {
    listData <- list(...)
    if (length(listData) && is(listData[[1L]], "MultiAssayExperiment"))
        multiassay <- listData[[1L]]
    else
        multiassay <- MultiAssayExperiment(...)
    multiassay <- .doMatching(multiassay)
    new("MatchedAssayExperiment", multiassay)
}

setAs("MultiAssayExperiment", "MatchedAssayExperiment", function(from) {
    from <- .doMatching(from)
    new("MatchedAssayExperiment", from)
})

