### ==============================================
### RangedRaggedAssay class
### ----------------------------------------------

#' An extension of the GRangesList class
#' 
#' @inheritParams GenomicRanges::GRangesList
#' @exportClass RangedRaggedAssay
.RangedRaggedAssay <- setClass("RangedRaggedAssay", contains = "GRangesList")

### - - - - - - - - - - - - - - - - - - - - - - - -
### Builder
###

#' Create a RangedRaggedAssay
#'
#' @param x A \code{list}, \code{GRanges} or \code{GRangesList} object
#' @return A \code{\linkS4class{RangedRaggedAssay}} class object
#' @export RangedRaggedAssay
RangedRaggedAssay <- function(x = GRangesList()){
  if(inherits(x, "GRanges")){
    x <- GRangesList(x)
  }
  .RangedRaggedAssay(x)
}
