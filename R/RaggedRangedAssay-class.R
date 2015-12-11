### ==============================================
### RaggedRangedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @inheritParams GenomicRanges::GRangesList
#' @exportClass RaggedRangedAssay
setClass("RaggedRangedAssay", contains = "GRangesList")

#' Generic Builder and Accessor Function
setGeneric("RaggedRangedAssay", function(x) standardGeneric("RaggedRangedAssay"))


### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert GRangesList to RaggedRangedAssay
#'
#' @param x A \code{GRangesList} object
#' @return A \code{\linkS4class{RaggedRangedAssay}} class object
#' @exportMethod RaggedRangedAssay
#' @describeIn RaggedRangedAssay Convert a GRangesList to a RaggedRangedAssay
setMethod("RaggedRangedAssay", "GRangesList", function(x) new("RaggedRangedAssay", x))
#' @describeIn RaggedRangedAssay Convert a GRanges to GRangesList to RaggedRangedAssay
setMethod("RaggedRangedAssay", "GRanges", function(x) new("RaggedRangedAssay", GRangesList(x)))
setMethod("RaggedRangedAssay", "ANY", function(x) x)
setMethod("RaggedRangedAssay", "list", function(x) lapply(x, RaggedRangedAssay))
