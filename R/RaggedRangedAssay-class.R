### ==============================================
### RaggedRangedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @inheritParams GenomicRanges::GRangesList
#' @exportClass RaggedRangedAssay
setClass("RaggedRangedAssay", contains = "GRangesList")

#' Generic Builder and Accessor Function
setGeneric("RaggedRangedAssay", function(...) standardGeneric("RaggedRangedAssay"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Convert GRangesList to RaggedRangedAssay
#'
#' @param ... A \code{GRangesList} object
#' @return A \code{\linkS4class{RaggedRangedAssay}} class object
#' @exportMethod RaggedRangedAssay
#' @describeIn RaggedRangedAssay Convert a GRangesList to a RaggedRangedAssay
setMethod("RaggedRangedAssay", "GRangesList", function(...){
  as(..., "RaggedRangedAssay")
})
setMethod("RaggedRangedAssay", "ANY", function(...){
  grl <- GRangesList(...)
  RaggedRangedAssay(grl)
})
