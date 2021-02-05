.MAE_RDS_BASENAME <- "mae.rds"
.EXPERIMENTS_H5_BASENAME <- "experiments.h5"

.load_HDF5Array_package <- function()
{
    if (!requireNamespace("HDF5Array", quietly=TRUE))
        stop("Please install the 'HDF5Array' package to use this function")
}

.isSingleString <- S4Vectors::isSingleString

.write_h5_experiments <- function(
    experiments, h5_path, chunkdim, level, as.sparse, verbose
) {
    HDF5Array:::.write_h5_assays(
        experiments, h5_path, chunkdim, level, as.sparse, verbose
    )
}

.serialize_HDF5SummarizedExperiment <- function(x, rds_path, verbose)
{
    assays(x) <- ExperimentList(.shorten_assay2h5_links(assays(x)))
    if (verbose)
        message("Serialize ", class(x), " object to ",
                ifelse(file.exists(rds_path), "existing ", ""),
                "RDS file:\n  ", rds_path)
    saveRDS(x, file=rds_path)
}

.write_HDF5MultiAssayExperiment <- function(x,
    rds_path=.MAE_RDS_BASENAME, h5_path=.EXPERIMENTS_H5_BASENAME,
    chunkdim=NULL, level=NULL, as.sparse=NA, verbose=FALSE
) {
    .load_HDF5Array_package()

    if (!is(x, "MultiAssayExperiment"))
        stop("<internal> 'x' must be a MultiAssayExperiment object")

    assays(x) <- .write_h5_assays(
        assays(x), h5_path, chunkdim, level, as.sparse, verbose
    )

    .seralize_HDF5MultiAssayExperiment(x, rds_path, verbose)
    invisible(x)
}

saveHDF5MultiAssayExperiment <-
    function(
        x, dir = tempdir(), prefix = NULL, replace = FALSE, chunkdim = NULL,
        level = NULL, as.sparse = NA, verbose = NA
    )
{
    if (is.null(prefix))
        prefix <- as.character(substitute(x))

    if (!requireNamespace("HDF5Array", quietly = TRUE))
        stop("Please install the 'HDF5Array' package to use this function")

    stopifnot(.isSingleString(dir), .isSingleString(prefix))

    if (!S4Vectors::isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")

    verbose <- DelayedArray:::normarg_verbose(verbose)

    if (!dir.exists(dir)) {
        .create_dir(dir)
    } else {
        .replace_dir(dir, replace)
    }
    rds_path <- file.path(dir, paste0(prefix, .SE_RDS_BASENAME))
    h5_path <- file.path(dir, paste0(prefix, .ASSAYS_H5_BASENAME))
    if (prefix != "")
        .check_and_delete_files(rds_path, h5_path, replace)

    .write_HDF5SummarizedExperiment(x, rds_path=rds_path,
                                       h5_path=h5_path,
                                       chunkdim=chunkdim, level=level,
                                       as.sparse=as.sparse,
                                       verbose=verbose)
}
