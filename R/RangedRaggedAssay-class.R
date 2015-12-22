### ==============================================
### RangedRaggedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @inheritParams GenomicRanges::GRangesList
#' @exportClass RangedRaggedAssay
setClass("RangedRaggedAssay", contains = "GRangesList")

#' Generic Builder and Accessor Function
setGeneric("RangedRaggedAssay", function(...) standardGeneric("RangedRaggedAssay"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert GRangesList to RangedRaggedAssay
#'
#' @param ... A \code{GRangesList} object
#' @return A \code{\linkS4class{RangedRaggedAssay}} class object
#' @exportMethod RangedRaggedAssay
#' @describeIn RangedRaggedAssay Convert a GRangesList to a RangedRaggedAssay
setMethod("RangedRaggedAssay", "GRangesList", function(...){
  as(..., "RangedRaggedAssay")
})
setMethod("RangedRaggedAssay", "ANY", function(...){
  grl <- GRangesList(...)
  RangedRaggedAssay(grl)
})
