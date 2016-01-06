### ==============================================
### MultiAssayView class
### ==============================================

#' An class used for staging a subset operation
#'
#' @slot subject A MultiAssayExperiment instance
#' @slot rowindex An \code{IntegerList} indexing rows of subject
#'     retained for subsequent work.
#' @slot colindex An \code{IntegerList} indexing columns of sujbect
#'     retained for subsequent work.
#' @exportClass MultiAssayView

.MultiAssayView <- setClass("MultiAssayView",
         representation(
             subject = "environment",
             rowindex="IntegerList",
             colindex="IntegerList")
         )

.idxlist <- function(x) {
    nms <- names(x)
    x <- unname(CharacterList(x))
    offset <- rep(cumsum(c(0L, lengths(x)[-length(x)])), lengths(x))
    idx <- seq_along(unlist(x, use.names=FALSE)) - offset
    setNames(relist(idx, x), nms)
}

MultiAssayView <- function(subject) {
    stopifnot(inherits(subject, "MultiAssayExperiment"))
    .MultiAssayView(subject=as.environment(list(subject=subject)),
                    rowindex=.idxlist(rownames(subject)),
                    colindex=.idxlist(colnames(subject)))
}

.subject <- function(x) {
    getElement(x, "subject")[["subject"]]
}

.rowindex <- function(x)
    getElement(x, "rowindex")

.colindex <- function(x)
    getElement(x, "colindex")

setMethod("rownames", "MultiAssayView", function(x) {
    CharacterList(rownames(.subject(x)))[.rowindex(x)]
})

setReplaceMethod("rownames", c("MultiAssayView", "ANY"),
    function(x, value)
{
    slot(x, "rownames") <-
        match(CharacterList(value), rownames(x))
    x
})

setMethod("colnames", "MultiAssayView", function(x) {
    CharacterList(colnames(.subject(x)))[.colindex(x)]
})

setReplaceMethod("colnames", c("MultiAssayView", "ANY"),
   function(x, value)
{
    slot(x, "colnames") <- CharacterList(x)
    x
})

.subset1 <- function(i, index, names)
    ## FIXME: below is +/- ok for typeof(i) == character only
    index[which(names %in% i)]

setMethod("[", c("MultiAssayView", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    rowindex <- .rowindex(x)
    if (!missing(i))
        rowindex <- .subset1(i, rowindex, rownames(x))
    colindex <- .colindex(x)
    if (!missing(j))
        colindex <- .subset1(j, colindex, colnames(x))
    initialize(x, rowindex=rowindex, colindex=colindex)
    ## FIXME: row/colnames should be updated to reflect subsets ?
})

materialize <- function(x)
{
    stopifnot(inherits(x, "MultiAssayView"))

    subject <- .subject(x)
    elist <- Elist(subject)
    elist <- mendoapply(function(e, i) e[i,,drop=FALSE], elist, rownames(x))
    elist <- mendoapply(function(e, j) e[,j, drop=FALSE], elist, colnames(x))
    initialize(subject, Elist=elist)
}

setMethod("show", "MultiAssayView", function(object) {
    rownames <- rownames(object)
    colnames <- colnames(object)
    cat("A", class(object), "object\n")
    cat("\nrownames():\n")
    print(rownames)
    cat("\ncolnames():\n")
    print(colnames)
})
