.MAE_RDS_BASENAME <- "mae.rds"
.EXPERIMENTS_H5_BASENAME <- "experiments.h5"

.load_HDF5Array_package <- function()
{
    if (!requireNamespace("HDF5Array", quietly=TRUE))
        stop("Please install the 'HDF5Array' package to use this function")
}

.isSingleString <- S4Vectors::isSingleString

.serialize_HDF5MultiAssayExperiment <- function(x, rds_path, verbose)
{
    exps <- HDF5Array::shorten_assay2h5_links(suppressWarnings(assays(x)))
    explist <- experiments(x)
    for (i in seq_along(exps)) {
        if (is(explist[[i]], "SummarizedExperiment")) {
            explist[[i]]@assays <- Assays(exps[[i]])
        } else if (is(explist[[i]], "HDF5Array")) {
            explist@listData[[i]] <- exps[[i]]
        }
    }
    x <- BiocBaseUtils::setSlots(x, ExperimentList = explist, check = FALSE)
    if (verbose)
        message("Serialize ", class(x), " object to ",
                ifelse(file.exists(rds_path), "existing ", ""),
                "RDS file:\n  ", rds_path)
    saveRDS(x, file=rds_path)
}

.write_h5_dimnames <- function(x, h5_path) {
    dimnamesList <- lapply(x, dimnames)
    h5_names <- sprintf("assay%03d", seq_along(x))
    invisible(
        Map(
            function(dl, al) {
                HDF5Array::h5writeDimnames(
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

    explist <- experiments(x)
    assaylist <- assays(x)
    exps <- HDF5Array::write_h5_assays(
        assaylist, h5_path, chunkdim, level, as.sparse, verbose
    )
    ## write Dimnames
    .write_h5_dimnames(assaylist, h5_path)
    ## avoid errors with arrays (withDimnames=FALSE)
    assays(explist, withDimnames = FALSE) <- exps
    x <- BiocBaseUtils::setSlots(
        x, ExperimentList = explist, check = FALSE
    )
    .serialize_HDF5MultiAssayExperiment(x, rds_path, verbose)
    invisible(x)
}

#' @rdname HDF5MultiAssayExperiment
#'
#' @title Save a MultiAssayExperiment class object to HDF5 and Rds files
#'
#' @description
#'   This function takes a `MultiAssayExperiment` object and uses the `assays`
#'   functionality to obtain data matrices out of the experiments. These are
#'   then saved into the `.h5` file format. This function relies heavily on
#'   the `HDF5Array` package whose installation is required before use.
#'   `saveHDF5MultiAssayExpeirment` preserves the classes contained in the
#'   \linkS4class{ExperimentList} with the exception of `matrix` which is
#'   converted to `HDF5Matrix`. Internal `SummarizedExperiment` assays are
#'   converted to HDF5-backed assays as in
#'   `HDF5Array::saveHDF5SummarizedExperiment`. `SummarizedExperiment`
#'   objects with multiple `i`-th assays will have the first assay take
#'   precedence and others assays will be dropped with a warning.
#'   If the first assay in a `SummarizedExperiment` contains an array,
#'   the array is preserved in the process of saving and loading the
#'   HDF5-backed `MultiAssayExperiment`.
#'
#' @inheritParams HDF5Array::saveHDF5SummarizedExperiment
#'
#' @param x A \linkS4class{MultiAssayExperiment} object or derivative
#'
#' @param dir The path (as a single string) to the directory where to save the
#'   HDF5-based \linkS4class{MultiAssayExperiment} object or to load it from.
#'
#'   When saving, the directory will be created if it doesn't already exist.
#'   If the directory already exists and no prefix is specified and
#'   `replace` is set to `TRUE`, then it's replaced with an
#'   empty directory.
#'
#' @param prefix An optional prefix to add to the names of the files created
#'   inside `dir`. This allows saving more than one object in the same
#'   directory. When the prefix is `NULL`, the name of the `x` input
#'   `MultiAssayExperiment` is used. To avoid the default setting
#'   use an empty character string i.e., `""`. An underscore (`_`) is
#'   appended to the prefix when provided; therefore, typical inputs should be
#'   words, e.g., "test".
#'
#' @param verbose Set to `TRUE` to make the function display progress.
#'
#'   In the case of `saveHDF5MultiAssayExperiment()`, `verbose`
#'   is set to `NA` by default, in which case verbosity is controlled
#'   by `DelayedArray.verbose.block.processing` option. Setting
#'   `verbose` to `TRUE` or `FALSE` overrides the option.
#'
#' @md
#'
#' @examples
#'
#' data("miniACC")
#'
#' testDir <- file.path(tempdir(), "test_mae")
#'
#' saveHDF5MultiAssayExperiment(
#'     miniACC, dir = testDir, verbose = TRUE, replace = TRUE
#' )
#'
#' ## inspect the files in the dir
#' list.files(testDir)
#'
#' loadHDF5MultiAssayExperiment(
#'     dir = testDir
#' )
#'
#' ## remove example files
#' unlink(testDir, recursive = TRUE)
#'
#' @export
saveHDF5MultiAssayExperiment <-
    function(
        x, dir = "h5_mae", prefix = NULL, replace = FALSE, chunkdim = NULL,
        level = NULL, as.sparse = NA, verbose = NA
    )
{
    if (is.null(prefix))
        prefix <- as.character(substitute(x))

    if (nzchar(prefix))
        prefix <- paste0(prefix, "_")

    .load_HDF5Array_package()

    stopifnot(.isSingleString(dir), .isSingleString(prefix))

    if (!BiocBaseUtils::isTRUEorFALSE(replace))
        stop("'replace' must be TRUE or FALSE")

    verbose <- .normarg_verbose(verbose)

    if (!dir.exists(dir)) {
        HDF5Array::create_dir(dir)
    } else {
        HDF5Array::replace_dir(dir, replace)
    }

    rds_path <- file.path(dir, paste0(prefix, .MAE_RDS_BASENAME))
    h5_path <- file.path(dir, paste0(prefix, .EXPERIMENTS_H5_BASENAME))

    if (prefix != "")
        HDF5Array::check_and_delete_files(rds_path, h5_path, replace)

    .write_HDF5MultiAssayExperiment(x,
        rds_path=rds_path, h5_path=h5_path, chunkdim=chunkdim, level=level,
        as.sparse=as.sparse, verbose=verbose
    )
}

#' @name HDF5MultiAssayExperiment
#'
#' @export
loadHDF5MultiAssayExperiment <- function(dir = "h5_mae", prefix = NULL)
{
    .load_HDF5Array_package()
    if (is.null(prefix)) {
        pattern <- paste0(.MAE_RDS_BASENAME, "|", .EXPERIMENTS_H5_BASENAME)
        prefix <- unique(gsub(pattern, "", dir(dir)))
        if (length(prefix) > 1L)
            stop("More than one object saved in 'dir', specify a 'prefix'")
    }
    if (nzchar(prefix) && !endsWith(prefix, "_"))
        prefix <- paste0(prefix, "_")

    stopifnot(.isSingleString(dir), .isSingleString(prefix))

    if (!dir.exists(dir)) {
        if (file.exists(dir))
            stop(wmsg("\"", dir, "\" is a file, not a directory"))
        stop(wmsg("directory \"", dir, "\" not found"))
    }

    rds_path <- file.path(dir, paste0(prefix, .MAE_RDS_BASENAME))
    ans <- try(.read_HDF5MultiAssayExperiment(rds_path), silent=TRUE)
    if (inherits(ans, "try-error"))
        HDF5Array::stop_if_bad_dir(dir, prefix)
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
    explist <- experiments(ans)

    for (i in seq_along(explist)) {
        if (is(explist[[i]], "SummarizedExperiment")) {
            explist[[i]]@assays <-
                HDF5Array::restore_absolute_assay2h5_links(
                    explist[[i]]@assays, dir
                )
        } else if (is(explist[[i]], "HDF5Array")) {
            explist <- S4Vectors::setListElement(explist, i,
                DelayedArray::modify_seeds(
                    explist[[i]],
                    function(x) {
                        what <- c("experiment ", i, " in the ",
                                  "MultiAssayExperiment object to load")
                        h5_path <- file.path(dir, x@filepath)
                        ## file_path_as_absolute() will fail if the file does
                        ## not exist.
                        if (!file.exists(h5_path))
                            stop(wmsg(what, " points to an HDF5 file ",
                                      "that does not exist: ", h5_path))
                        x@filepath <- tools::file_path_as_absolute(h5_path)
                        ## Check that 'x' points to an HDF5 dataset that is
                        ## accessible and "as expected".
                        msg <- HDF5Array::validate_HDF5ArraySeed_dataset_geometry(x, what)
                        if (!isTRUE(msg))
                            stop(wmsg(msg))
                        x
                    }
                )
            )
        }
    }
    ans <- BiocBaseUtils::setSlots(
        ans, ExperimentList = explist, check = FALSE
    )
    ans <- updateObject(ans, check=FALSE)
    if (!is(ans, "MultiAssayExperiment"))
        stop(wmsg("the object serialized in \"", rds_path, "\" is not ",
                  "a MultiAssayExperiment object or derivative"))
    ans
}
