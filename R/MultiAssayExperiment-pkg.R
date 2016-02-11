#' MultiAssayExperiment: Build an integrative multiassay container
#' 
#' MultiAssayExperiment allows the manipulation of related multiassay 
#' datasets with partially overlapping samples, associated metadata at 
#' the level of an entire study, and at the level of the "biological unit".
#' The biological unit may be a patient, plant, yeast strain, etc.
#'
#' The package hierarch of information: 
#' \itemize{
#' \item study
#' \item experiments
#' \item samples
#' }
#' 
#' @docType package
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom Biobase pData
#' @importFrom IRanges CharacterList
"_PACKAGE"