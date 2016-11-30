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
#' @param ... Additional arguments, see details for more information.
#'
#' @details
#' The \linkS4class{RangedRaggedAssay} class represents genomic ranges in
#' matrix shape based on the selected inner metadata column ("mcol"). To
#' accomplish this, the \code{mcolname} argument can be indicated with a
#' string. The \code{background} argument can be used to specify a background
#' value for the resulting matrix (default NA). This usually indicates non-matching
#' values in the matrix (e.g., 2 for diploid genomes). Users are also able
#' to provide a \link{GRanges} class object for specifying ranges of
#' interest in the resulting matrix using the \code{ranges} argument.
#' The \code{make.names} argument is a logical value (default FALSE) that
#' allows the user to indicate the automatic creation of automatic names
#' either from the \code{GRanges} or the \link{RangedRaggedAssay} object in
#' character format (i.e., "chr1:2-3:+"). The user can also include a
#' \code{\strong{type}} argument for indicating the type of overlap check
#' requested (default "any").
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
#' @exportMethod assay
setMethod("assay", c("RangedRaggedAssay", "missing"),
          function(x, i, ...) {
              args <- list(...)
              if (!all(GenomicRanges::isDisjoint(x)))
                  stop("only disjoint ranges supported")

              if (is.null(args$mcolname))
                  args$mcolname <- "score"
              if (!is.numeric(mcols(x[[1L]])[[args$mcolname]]))
                  stop("metadata column is not numeric")
              if (is.null(args$background))
                  args$background <- NA
              if (is.null(args$make.names))
                  args$make.names <- FALSE
              if (is.null(args$type))
                  args$type <- "any"

              if (!is.null(args$ranges)) {
                  ranges <- args$ranges
                  if (!inherits(ranges, "GRanges"))
                      stop("ranges must be a GRanges object")
                  if (args$make.names || is.null(names(ranges))) {
                      rowNames <- as.character(ranges)
                  } else {
                      rowNames <- names(ranges)
                  }
              } else {
                  rowNames <- rownames(x)
                  if (args$make.names) {
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
                                                              type = args$type),
                                         as.numeric(mcols(
                                             obj[[j]])[[args$mcolname]]
                                         ),
                                         args$background)
                                     return(MValues)
                                 }, obj = x))
              colnames(newMatrix) <- colnames(x)
              rownames(newMatrix) <- rowNames
              return(newMatrix)
          })

#' @describeIn ExperimentList Get the assay data for the default ANY class
setMethod("assay", c("ANY", "missing"), function(x, i) {
    if (inherits(x, "ExpressionSet"))
        return(Biobase::exprs(x))
    I(x)
})

#' @describeIn ExperimentList Get the assay data from each element in the
#' \link{ExperimentList}
#' @param i missing argument
#' @aliases assay,ExperimentList,missing-method
setMethod("assay", c("ExperimentList", "missing"), function(x, i) {
    lapply(x, FUN = function(y) {assay(y)})
})

#' @describeIn MultiAssayExperiment Get the assay data for a
#' \link{MultiAssayExperiment} as a \code{list}
#' @aliases assay,MultiAssayExperiment,missing-method
setMethod("assay", c("MultiAssayExperiment", "missing"), function(x, i) {
    assay(experiments(x))
})
