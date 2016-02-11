### ==============================================
### MultiAssayView class
### ==============================================

#' A class used for staging a subset operation
#'
#' Use a \code{MultiAssayView} class to initialize a subsetting operation
#' of a \code{\linkS4class{MultiAssayExperiment}} class object
#'
#' @slot subject A \code{\linkS4class{MultiAssayExperiment}} instance
#' @slot rowindex An \code{IntegerList} indexing rows of subject
#'     retained for subsequent operations
#' @slot colindex An \code{IntegerList} indexing columns of sujbect
#'     retained for subsequent operations
#' @slot assayindex An \code{integer} indexing assays of subject
#'     retained for subsequent operations
#' @return A \code{MultiAssayView} object
#' @exportClass MultiAssayView
#' @aliases MultiAssayView
.MultiAssayView <- setClass("MultiAssayView",
         representation(
             subject = "environment",
             rowindex = "IntegerList",
             colindex = "IntegerList", 
             assayindex = "integer"
         ))

.idxlist <- function(x) {
    nms <- names(x)
    x <- unname(CharacterList(x))
    offset <- rep(cumsum(c(0L, lengths(x)[-length(x)])), lengths(x))
    idx <- seq_along(unlist(x, use.names=FALSE)) - offset
    stats::setNames(relist(idx, x), nms)
}

#' @export MultiAssayView
MultiAssayView <- function(subject) {
    stopifnot(inherits(subject, "MultiAssayExperiment"))
    .MultiAssayView(subject=as.environment(list(subject=subject)),
                    rowindex=.idxlist(rownames(subject)),
                    colindex=.idxlist(colnames(subject)),
                    assayindex = seq_along(subject))
}

.subject <- function(x) {
    getElement(x, "subject")[["subject"]]
}

.rowindex <- function(x)
    getElement(x, "rowindex")

.colindex <- function(x)
    getElement(x, "colindex")

.assayindex <- function(x)
    getElement(x, "assayindex")

#' @describeIn MultiAssayView Get a CharacterList of rownames
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

#' @describeIn MultiAssayView Get a CharacterList of colnames
setMethod("colnames", "MultiAssayView", function(x) {
    CharacterList(colnames(.subject(x)))[.colindex(x)]
})

setReplaceMethod("colnames", c("MultiAssayView", "ANY"),
   function(x, value)
{
    slot(x, "colnames") <- CharacterList(x)
    x
})

#' @describeIn MultiAssayView Get a \code{character} vector of experiment names
setMethod("names", "MultiAssayView", function(x)
  names(getElement(x, "subject")[["subject"]])
)

#' @describeIn MultiAssayView Get the number of assays in the
#' \code{MultiAssayExperiment} instance
setMethod("length", "MultiAssayView", function(x)
  length(getElement(x, "subject")[["subject"]])
)

.subset1 <- function(i, index, names)
    ## FIXME: below is +/- ok for typeof(i) == character only
    index[which(names %in% i)]

#' @describeIn MultiAssayView Subset MultiAssayView dimensions with a 
#' \code{character} vector
#' @param x A \code{MultiAssayView} class object
#' @param i A \code{character} vector for subsetting rownames
#' @param j A \code{character} vector for subsetting colnames
#' @param k A \code{character} vector for subsetting assays
#' @param ... Additional parameters 
#' @param drop logical (default TRUE) whether to drop empty assay elements
setMethod("[", c("MultiAssayView", "ANY", "ANY", "ANY"),
    function(x, i, j, k, ..., drop=TRUE)
{
    rowindex <- .rowindex(x)
    if (!missing(i))
        rowindex <- .subset1(i, rowindex, rownames(x))
    colindex <- .colindex(x)
    if (!missing(j))
        colindex <- .subset1(j, colindex, colnames(x))
    assayindex <- .assayindex(x)
    if (!missing(k))
        assayindex <- assayindex[k]
    initialize(x, rowindex=rowindex, colindex=colindex, assayindex = assayindex)
    ## FIXME: row/colnames should be updated to reflect subsets ?
})

materialize <- function(x)
{
    stopifnot(inherits(x, "MultiAssayView"))

    subject <- .subject(x)
    elist <- Elist(subject)
    elist <- elist[names(x)]
    elist <- mendoapply(function(e, i) e[i,,drop=FALSE], elist, rownames(x))
    elist <- mendoapply(function(e, j) e[,j, drop=FALSE], elist, colnames(x))
    initialize(subject, Elist=elist)
}

#' @describeIn MultiAssayView Show method for \code{MultiAssayView}
#' @param object A \code{MultiAssayView} class object
setMethod("show", "MultiAssayView", function(object) {
    rownames <- rownames(object)
    colnames <- colnames(object)
    cat("A", class(object), "object\n")
    cat("\nrownames():\n")
    print(rownames)
    cat("\ncolnames():\n")
    print(colnames)
})
