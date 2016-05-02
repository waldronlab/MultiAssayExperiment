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
#' @param ranges A \link{GRanges} class identifying the ranges of interest
#' @param background A single value for the non-matching background values in
#' the matrix (e.g., 2 for diploid genomes)
#' @param use.names logical (default FALSE) whether to use names in the
#' given ranges argument or in the provided 'x' argument
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
setMethod("assay", "RangedRaggedAssay",
          function(x, ranges = NULL, background = NA, use.names = FALSE) {
  if (!all(GenomicRanges::isDisjoint(x))) {
    stop("Matrix can only be created for disjoint ranges")
  }
  if (.uniqueSortIdentical(names(x), names(metadata(x))) &&
      !identical(names(metadata(x)), names(x))) {
    metadata(x) <- metadata(x)[match(names(x), names(metadata(x)))]
  }
  if (!is.null(ranges)) {
    if (!inherits(ranges, "GRanges")) {
      stop("ranges must be a GRanges object")
    }
    if (use.names) {
      rowNames <- names(ranges)
    } else {
      rowNames <- as.character(ranges)
    }
  } else {
    rangeNames <- unique(as.character(unlist(x, use.names = FALSE)))
    ranges <- as(rangeNames, "GRanges")
    if (use.names) {
      featNames <- names(unlist(x, use.names = FALSE))
      if (length(unique(featNames)) == length(rangeNames)) {
        rowNames <- featNames
      } else {
        stop("feature names not unique accross ranges")
      }
    }
    rowNames <- rangeNames
  }
  newMatrix <- do.call(cbind, lapply(seq_along(x), function(i, obj) {
    MValues <- ifelse(
      IRanges::overlapsAny(ranges, obj[[i]], type = "equal"), 
      as.numeric(metadata(obj)[[i]]$score),
      background
    )
    return(MValues)
  }, obj = x))
  colnames(newMatrix) <- names(x)
  rownames(newMatrix) <- rowNames
  return(newMatrix)
})

#' @describeIn Elist Get the assay data for the default ANY class
setMethod("assay", "ANY", function(x) {
  I(x)
})

#' @describeIn Elist Get the assay data from each element in the \link{Elist}
setMethod("assay", "Elist", function(x) {
  lapply(x, FUN = function(y) {assay(y)})
})

#' @describeIn MultiAssayExperiment Get the assay data for a
#' \link{MultiAssayExperiment} as a \code{list}
setMethod("assay", "MultiAssayExperiment", function(x) {
  assay(Elist(x))
})
