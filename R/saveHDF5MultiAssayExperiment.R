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
    exps <- as(.shorten_assay2h5_links(assays(x)), "ExperimentList")
    x <- BiocGenerics:::replaceSlots(x, ExperimentList = exps, check = FALSE)
    if (verbose)
        message("Serialize ", class(x), " object to ",
                ifelse(file.exists(rds_path), "existing ", ""),
                "RDS file:\n  ", rds_path)
    saveRDS(x, file=rds_path)
}

.write_h5_dimnames <- function(x, h5_path) {
    dimnamesList <- lapply(experiments(x), dimnames)
    h5_names <- sprintf("assay%03d", seq_along(x))
    invisible(
        Map(
            function(dl, al) {
                HDF5Array:::h5writeDimnames(
                    dimnames = dl, filepath = h5_path, name = al,
                    group = NA, h5dimnames = NULL
                )
            },
            dl = dimnamesList, al = h5_names
        )
    )
}

.write_HDF5MultiAssayExperiment <- function(x,
    rds_path=.MAE_RDS_BASENAME, h5_path=.EXPERIMENTS_H5_BASENAME,
    chunkdim=NULL, level=NULL, as.sparse=NA, verbose=FALSE
) {
    .load_HDF5Array_package()

    if (!is(x, "MultiAssayExperiment"))
        stop("<internal> 'x' must be a MultiAssayExperiment object")

    exps <- HDF5Array:::.write_h5_assays(
        assays(x), h5_path, chunkdim, level, as.sparse, verbose
    )
    ## write Dimnames
    .write_h5_dimnames(x, h5_path)

    experiments(x) <- exps
    .serialize_HDF5MultiAssayExperiment(x, rds_path, verbose)
    invisible(x)
}

#' Save a MultiAssayExperiment class object to HDF5 and Rds files
#'
#' This function takes a `MultiAssayExperiment` object and uses the `assays`
#' functionality to create data matrices out of ti
#'
#' @rd
#'
#' @examples
#'
#' testDir <- tempdir()
#' saveHDF5MultiAssayExperiment(
#'     miniACC, dir = file.path(tempdir(), "test_mae"), replace = TRUE
#' )
#'
#'
#' @export
saveHDF5MultiAssayExperiment <-
    function(
        x, dir = tempdir(), prefix = NULL, replace = FALSE, chunkdim = NULL,
        level = NULL, as.sparse = NA, verbose = NA
    )
{
    if (is.null(prefix))
        prefix <- paste0(as.character(substitute(x)), "_")

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

loadHDF5MultiAssayExperiment <- function(dir = "h5_mae", prefix = NULL)
{
    .load_HDF5Array_package()

    if (is.null(prefix))
        prefix <- unique(
            vapply(strsplit(dir(dir), "_"), '[', character(1L), 1L)
        )

    stopifnot(.isSingleString(dir), .isSingleString(prefix))

    if (!dir.exists(dir)) {
        if (file.exists(dir))
            stop(wmsg("\"", dir, "\" is a file, not a directory"))
        stop(wmsg("directory \"", dir, "\" not found"))
    }

    rds_path <- file.path(dir, paste0(prefix, .MAE_RDS_BASENAME))
    ans <- try(.read_HDF5MultiAssayExperiment(rds_path), silent=TRUE)
    if (inherits(ans, "try-error"))
        HDF5Array:::.stop_if_bad_dir(dir, prefix)
    ans
}

.read_HDF5MultiAssayExperiment <- function(rds_path)
{
    if (!file.exists(rds_path))
        stop(wmsg("file not found: ", rds_path))
    if (dir.exists(rds_path))
        stop(wmsg("'", rds_path, "' is a directory, not a file"))

    dir <- dirname(rds_path)
    ans <- readRDS(rds_path)
    ans@ExperimentList <- HDF5Array:::.restore_absolute_assay2h5_links(
        ans@ExperimentList, dir
    )
    ans <- updateObject(ans, check=FALSE)
    if (!is(ans, "MultiAssayExperiment"))
        stop(wmsg("the object serialized in \"", rds_path, "\" is not ",
                  "a MultiAssayExperiment object or derivative"))
    ans
}
