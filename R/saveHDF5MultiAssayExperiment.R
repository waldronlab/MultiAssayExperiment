.MAE_RDS_BASENAME <- "mae.rds"
.EXPERIMENTS_H5_BASENAME <- "experiments.h5"

.load_HDF5Array_package <- function()
{
    if (!requireNamespace("HDF5Array", quietly=TRUE))
        stop("Please install the 'HDF5Array' package to use this function")
}

.isSingleString <- S4Vectors::isSingleString
.shorten_assay2h5_links <- HDF5Array:::.shorten_assay2h5_links

.serialize_HDF5MultiAssayExperiment <- function(x, rds_path, verbose)
{
    experiments(x) <- ExperimentList(.shorten_assay2h5_links(assays(x)))
    if (verbose)
        message("Serialize ", class(x), " object to ",
                ifelse(file.exists(rds_path), "existing ", ""),
                "RDS file:\n  ", rds_path)
    saveRDS(x, file=rds_path)
}

.write_h5_dimnames <- function(x, h5_path) {
    dimnamesList <- lapply(experiments(x), dimnames)
    h5_names <- sprintf("assay%03d", seq_along(x))
    Map(
        function(x, y) {
            HDF5Array:::h5writeDimnames(
                dimnames = x, filepath = h5_path, name = y,
                group = NA, h5dimnames = NULL
            )
        },
        x = dimnamesList, y = h5_names
    )
}

.write_HDF5MultiAssayExperiment <- function(x,
    rds_path=.MAE_RDS_BASENAME, h5_path=.EXPERIMENTS_H5_BASENAME,
    chunkdim=NULL, level=NULL, as.sparse=NA, verbose=FALSE
) {
    .load_HDF5Array_package()

    if (!is(x, "MultiAssayExperiment"))
        stop("<internal> 'x' must be a MultiAssayExperiment object")

    experiments(x) <- HDF5Array:::.write_h5_assays(
        assays(x), h5_path, chunkdim, level, as.sparse, verbose
    )
    ## write Dimnames
    .write_h5_dimnames(x, h5_path)

    .serialize_HDF5MultiAssayExperiment(x, rds_path, verbose)
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

    .load_HDF5Array_package()

    stopifnot(.isSingleString(dir), .isSingleString(prefix))

    if (!S4Vectors::isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")

    verbose <- DelayedArray:::normarg_verbose(verbose)

    if (!dir.exists(dir)) {
        HDF5Array:::.create_dir(dir)
    } else {
        HDF5Array:::.replace_dir(dir, replace)
    }

    rds_path <- file.path(dir, paste0(prefix, .MAE_RDS_BASENAME))
    h5_path <- file.path(dir, paste0(prefix, .EXPERIMENTS_H5_BASENAME))

    if (prefix != "")
        HDF5Array:::.check_and_delete_files(rds_path, h5_path, replace)

    .write_HDF5MultiAssayExperiment(x,
        rds_path=rds_path, h5_path=h5_path, chunkdim=chunkdim, level=level,
        as.sparse=as.sparse, verbose=verbose
    )
}
