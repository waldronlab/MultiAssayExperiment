#' MultiAssayExperiment: Build an integrative multi-assay container
#'
#' MultiAssayExperiment allows the manipulation of related multiassay
#' datasets with partially overlapping samples, associated metadata at
#' the level of an entire study, and at the level of the "biological unit".
#' The biological unit may be a patient, plant, yeast strain, etc.
#'
#' The package hierarchy of information:
#' \itemize{
#' \item study
#' \item experiments
#' \item samples
#' }
#'
#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom Biobase pData
#' @importFrom IRanges CharacterList
#' @aliases NULL
"_PACKAGE"

