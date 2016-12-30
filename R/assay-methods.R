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
#' @param i Argument from generic (default 1L)
#' @param mcolname A single string indicating the metadata column to use for
#' the assay conversion
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
#'
#' @examples
#' example("RangedRaggedAssay")
#'
#' ## Add some phony metadata to the RangedRaggedAssay
#' metadata(myRRA) <- list(snparrray1 = DataFrame(score = 1),
#' snparray2 = DataFrame(score = 1),
#' snparray3 = DataFrame(score = 3))
#'
#' assay(myRRA, background = 2)
#'
#' @return A \code{matrix} of values from the score column of the metadata.
#' @seealso \link{overlapsAny}
#' @exportMethod assay
setMethod("assay", c("RangedRaggedAssay", "missing"),
          function(x, i, mcolname = "score", background = NA,
                   make.names = FALSE, ranges = NULL, type = "any"){
              if (!all(GenomicRanges::isDisjoint(x)))
                  stop("only disjoint ranges supported")

              if (!is.numeric(mcols(x[[1L]])[[mcolname]]))
                  stop("metadata column is not numeric")

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
                  ranges <- GenomicRanges::GRanges(unlist(x, use.names = FALSE))
              }
              newMatrix <-
                  do.call(cbind,
                          lapply(seq_along(x),
                                 function(j, obj) {
                                     MValues <- ifelse(
                                         IRanges::overlapsAny(ranges, obj[[j]],
                                                              type = type),
                                         as.numeric(mcols(
                                             obj[[j]])[[mcolname]]
                                         ),
                                         background)
                                     return(MValues)
                                 }, obj = x))
              colnames(newMatrix) <- colnames(x)
              rownames(newMatrix) <- rowNames
              return(newMatrix)
          })

#' @describeIn ExperimentList Get the assay data for the default ANY class
setMethod("assay", c("ANY", "missing"), function(x, i) {
    if (is(x, "ExpressionSet"))
        return(Biobase::exprs(x))
    I(x)
})

#' @describeIn ExperimentList Get the assay data from each element in the
#' \link{ExperimentList}
#' @param i missing argument
#' @aliases assay,ExperimentList,missing-method
setMethod("assay", c("ExperimentList", "missing"), function(x, i) {
    lapply(x, FUN = function(y) { assay(y) })
})

#' @describeIn MultiAssayExperiment Get the assay data for a
#' \link{MultiAssayExperiment} as a \code{list}
#' @aliases assay,MultiAssayExperiment,missing-method
setMethod("assay", c("MultiAssayExperiment", "missing"), function(x, i) {
    assay(experiments(x))
})
