#' Subset MultiAssayExperiment object
#'
#' The \code{subset} function allows for searching through colnames and
#' rownames to obtain the requested subset of data or experiments (assays).
#'
#' @param x A \code{\linkS4class{MultiAssayExperiment}} object
#' @param indicator A \code{logical} or \code{character} vector or
#' \code{GRanges} class object to use for subsetting
#' @param method A \code{character} vector of length one designating to subset
#' either by columns, rows, or assays
#' @param drop logical (default FALSE) whether to coerce lowest possible
#' dimension after subsetting
#' @param ... Additional arguments to pass to
#' \code{\link[IRanges]{subsetByOverlaps}} when subsetting
#' by rownames
#' @return A subsetted \link{MultiAssayExperiment} class object
#' @exportMethod subset
#' @seealso subsetByAssay, subsetByColumn, subsetByRow
setMethod("subset", "MultiAssayExperiment",
          function(x, indicator, method = NULL, drop = TRUE, ...) {
            if (inherits(indicator, "GRanges")) {
              method <- "rows"
            } else if (is.null(method)) {
              stop("Indicate a subset method")
            } else {
              method <- match.arg(method, c("columns", "rows", "assays"))
            }
            if (method == "columns") {
              MultiAssay <- subsetByColumn(x = x,
                                           y = indicator)
            } else if (method == "rows") {
              MultiAssay <- subsetByRow(x = x,
                                        y = indicator, ...)
            } else if (method == "assays") {
              MultiAssay <- subsetByAssay(x = x,
                                          y = indicator)
            }
            if (drop) {
              isEmptyAssay <- vapply(Elist(MultiAssay), FUN = .isEmpty, FUN.VALUE = logical(1L))
              if (all(isEmptyAssay)) {
                warning("no data in assays")
                Elist(MultiAssay) <- Elist()
              } else if (any(isEmptyAssay)) {
                keeps <- names(isEmptyAssay)[sapply(isEmptyAssay, function(z) !isTRUE(z))]
                MultiAssay <- MultiAssay[, , keeps, drop = FALSE]
              }
            }
            return(MultiAssay)
          })
