#' Check expression of a given feature against clinical variable
#'
#' Function that outputs a \link{DataFrame} with participant ID, sample ID,
#' the select colData column, the expression values for select rownames,
#' and the center values for each gene by cluster.
#'
#' @param MultiAssayExperiment A MultiAssayExperiment object
#' @param colDataCols Select columns from the MultiAssayExperiment colData
#' DataFrame
#' @param rownames Features to be used for clustering
#' (e.g., a set of gene names)
#' @param experiments A \code{character} vector indicating assays of interest
#' in the \code{ExperimentList}
#' @param seed A single integer value passed to \link{set.seed} (default NULL)
#'
#' @return A DataFrame with appended cluster and center values
#' @examples
#' example(MultiAssayExperiment)
#' clusterOn(myMultiAssayExperiment, colDataCols = "sex",
#'     rownames = c("XIST", "RPS4Y1", "KDM5D"),
#'     experiments = "RNASeqGene", seed = 42L)
#'
#' @export clusterOn
#' @importFrom stats kmeans
clusterOn <- function(MultiAssayExperiment, colDataCols, rownames,
                      experiments, seed = NULL) {
    MultiAssayExperiment <- MultiAssayExperiment[rownames, , experiments]
    longMulti <- rearrange(MultiAssayExperiment, colDataCols = colDataCols)

    wideMulti <- tidyr::spread(
        as.data.frame(longMulti)[, -(which(names(longMulti)=="assay"))],
        key = "rowname", value = "value")

    expSubset <- wideMulti[, rownames]
    scaledExp <- apply(expSubset, 2, scale)
    if (!is.null(seed))
        set.seed(seed)
    kma <- kmeans(scaledExp, centers = 2)
    kcenters <- kma[["centers"]][kma[["cluster"]], ]
    rownames(kcenters) <- NULL
    colnames(kcenters) <- paste0("centers.", colnames(kcenters))
    kcenters <- cbind(kcluster = kma[["cluster"]], kcenters)

    sexcomp <- S4Vectors::DataFrame(wideMulti, kcenters)
    return(sexcomp)
}
