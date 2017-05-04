#' Create a Matrix of score values using a GRanges or own ranges
#'
#' This function can take a GRanges argument and use each range to check for
#' overlaps with any of the current ranges in the first argument and return a
#' score value from the corresponding metadata. This function will only operate
#' on fully disjoint ranges (see \code{isDisjoint} for details). It can only
#' work if metadata is present and there is a "score" column in the metadata.
#' Please see example on how to add metadata to a
#' \linkS4class{RangedRaggedAssay} or \link{GRangesList} class. This function
#' uses the \link{overlapsAny} function from the \code{GenomicRanges} package.
#'
#' @param x A \linkS4class{RangedRaggedAssay} or \link{GRangesList} class
#' @param i Argument set to missing (not used)
#' @param mcolname A single string indicating the metadata column to use for
#' the values in the resulting assay matrix
#' @param background A default background value for the resulting assay matrix
#' (default NA). This works for non-matching sample and range pairs in the data
#' and will be imputed in the matrix (e.g., 2 for diploid genomes)
#' @param make.names logical (default FALSE) whether to create character format
#' ranges for the rows of the matrix (either from the \code{ranges} argument
#' or from the \code{RangedRaggedAssay} itself). Example character format:
#' "chr1:2-3:+"
#' @param ranges An optional \link{GRanges} object for comparing accross all
#' sample ranges and for superseding the rows for the resulting matrix
#' (default NULL)
#' @param type The type argument from \link{overlapsAny}
#' @param ... Unused argument
#'
#' @return A \code{matrix} of values from the score column of the metadata.
#' @seealso \link{overlapsAny}
#' @exportMethod assay
setMethod("assay", c("RangedRaggedAssay", "missing"),
          function(x, i, mcolname = "score", background = NA,
                   make.names = FALSE, ranges = NULL, type = "any", ...) {
              .Defunct("RaggedExperiment")
              if (!all(GenomicRanges::isDisjoint(x)))
                  stop("only disjoint ranges supported")
              if (!is.null(ranges)) {
                  if (!is(ranges, "GRanges"))
                      stop("ranges must be a GRanges object")
                  if (make.names || is.null(names(ranges))) {
                      rowNames <- as.character(ranges)
                  } else {
                      rowNames <- names(ranges)
                  }
              } else {
                  rowNames <- rownames(x)
                  if (make.names) {
                      rangeNames <- unique(as.character(
                          unlist(x, use.names = FALSE)))
                      if (length(unique(rowNames)) != length(rangeNames))
                          stop("feature names not unique accross ranges")
                      rowNames <- rangeNames
                  }
                  ranges <- GenomicRanges::GRanges(unlist(x,
                                                          use.names = FALSE))
              }
              newMatrix <-
                  do.call(cbind,
                          lapply(seq_along(x),
                                 function(j, obj) {
                                     MValues <- ifelse(
                                         IRanges::overlapsAny(ranges, obj[[j]],
                                                              type = type),
                                         mcols(obj[[j]])[[mcolname]],
                                         background)
                                     return(MValues)
                                 }, obj = x))
              colnames(newMatrix) <- colnames(x)
              rownames(newMatrix) <- rowNames
              return(newMatrix)
          })

#' @describeIn ExperimentList Get the assay data for the default ANY class
setMethod("assay", c("ANY", "missing"), function(x, i, ...) {
    if (is(x, "ExpressionSet"))
        return(Biobase::exprs(x))
    return(x)
})


#' @describeIn ExperimentList Get the assay data from each element in the
#' \link{ExperimentList}
#' @param withDimnames logical (default TRUE) whether to return dimension names
#' @aliases assay,ExperimentList,missing-method
setMethod("assays", "ExperimentList", function(x, ..., withDimnames = TRUE) {
    as(IRanges::endoapply(x, FUN = function(y) assay(y, ...)), "SimpleList")
})

#' @describeIn ExperimentList Convenience function for the assay of the first
#' element
#' @param i A scalar \code{character} or \code{integer} index
setMethod("assay", c("ExperimentList", "missing"), function(x, i, ...) {
    if (!length(x))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(<", class(x), ">) is 0'")
    assay(x[[1L]], ...)
})

#' @describeIn ExperimentList Obtain the specified assay from ExperimentList
#' with a \code{numeric} index
setMethod("assay", c("ExperimentList", "numeric"), function(x, i, ...) {
    tryCatch({
        assay(x[[i]], ...)
    }, error = function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

#' @describeIn ExperimentList Get the specified assay from ExperimentList with
#' a \code{character} index
setMethod("assay", c("ExperimentList", "character"), function(x, i, ...) {
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assay(x[[i]], ...)
    }, error = function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(<", class(x), ">)")
    res
})

#' @describeIn MultiAssayExperiment Obtain a \code{\link{SimpleList}} of assay
#' data for all available experiments in the object
#' @param withDimnames logical (default TRUE) whether to return dimension names
#' included in the object
#' @exportMethod assays
setMethod("assays", "MultiAssayExperiment", function(x, ..., withDimnames = TRUE) {
    assays(experiments(x), ..., withDimnames = withDimnames)
})

#' @describeIn MultiAssayExperiment Convenience function for extracting the
#' assay of the first element in the \code{ExperimentList}
setMethod("assay", c("MultiAssayExperiment", "missing"), function(x, i, ...) {
    assay(experiments(x), ...)
})

#' @describeIn MultiAssayExperiment Obtain the specified assay from the
#' MultiAssayExperiment with a \code{numeric} index
setMethod("assay", c("MultiAssayExperiment", "numeric"), function(x, i, ...) {
    assay(experiments(x), i = i, ...)
})

#' @describeIn MultiAssayExperiment Get the specified assay from the
#' MultiAssayExperiment with a \code{character} index
setMethod("assay", c("MultiAssayExperiment", "character"), function(x, i, ...) {
    assay(experiments(x), i = i, ...)
})
