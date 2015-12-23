### ==============================================
### RangedRaggedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @inheritParams GenomicRanges::GRangesList
#' @exportClass RangedRaggedAssay
setClass("RangedRaggedAssay", contains = "GRangesList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

.RangedRaggedAssay <- setClass("RangedRaggedAssay", contains = "GRangesList")

#' Create a RangedRaggedAssay
#'
#' @param x A \code{list}, \code{GRanges} or \code{GRangesList} object
#' @return A \code{\linkS4class{RangedRaggedAssay}} class object
#' @export RangedRaggedAssay
RangedRaggedAssay <- function(x = list()){
  .RangedRaggedAssay(GRangesList(x))
}
