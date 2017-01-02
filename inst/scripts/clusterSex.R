# Check sex expression against clinical sex
#
# Function that outputs a data.frame with participant ID, sex,
# expression value, participant's sex from pData,
# the expression value for select genes, and the center value for each gene
# by cluster (matrix of cluster centers).
#
# @params MultiAssayExperiment A MultiAssayExperiment object
# @params pDataCols Select columns from the MultiAssayExperiment pData
# @params  ExperimentList.rows Genes to be used
# @params ExperimentList.exp A ExperimentList class object of experiment data
# @export clusterSex
clusterSex <- function(MultiAssayExperiment, pDataCols, ExperimentList.rows,
                      ExperimentList.exp, seed) {
    #Subset the MAE by rows(genes) and experiments
    MAE_s <- MAE[ExperimentList.rows, , ExperimentList.exp]
    #Use gather to get a dataframe with the metadata and experiment data
    meta_MAE <- gather(MAE_s, pDataCols = pDataCols)

    ##LONG TO WIDE
    meta_MAEw <- reshape(as.data.frame(meta_MAE), timevar = "rowname",
                         idvar = c("primary", "colname", pDataCols),
                         drop ="assay", direction = "wide")

    #Normalization of values for each experimentList.
    subset <- meta_MAEw[, sapply(meta_MAEw, is.numeric)]
    temp <- matrix(NA, nrow = nrow(meta_MAEw),
                  ncol = length(ExperimentList.rows)) # temp as a matrix
    #browser()
    for (i in seq_along(subset)) {
        temp[, i] <- scale(subset[, i])
    }

    meta_MAEw_n <- as.data.frame(temp)
    colnames(meta_MAEw_n) <- ExperimentList.rows
    meta_MAEw_n$PatientID <- meta_MAEw$primary
    meta_MAEw_n$sampleID <- meta_MAEw$colname
    meta_MAEw_n$Mgender <- meta_MAEw$gender

    set.seed(seed)
    kma <- kmeans(meta_MAEw_n[, ExperimentList.rows], centers = 2)
    meta_MAEw_n$Cluster <- kma$cluster
    centers_a <- kma$centers
    Kcenters <- centers_a[meta_MAEw_n[["Cluster"]], ]
    colnames(Kcenters) <- ExperimentList.rows

    genderscomp <- cbind(meta_MAEw_n[, c("PatientID", "Mgender", "Cluster")],
                         Kcenters)

    return(genderscomp)
}
