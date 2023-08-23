#' @include MultiAssayExperiment-subset.R
NULL

.warnnomatch <- function(type, values) {
    warning("'", type, "': ",
        paste(BiocBaseUtils::selectSome(values), collapse = ", "),
            "\n  could not be matched")
}

.sampleMapFromData <- function(colData, experiments) {
    samps <- colnames(experiments)
    assay <- factor(rep(names(samps), lengths(samps)), levels=names(samps))
    colname <- unlist(samps, use.names=FALSE)
    matches <- match(colname, rownames(colData))
    if (length(matches) && all(is.na(matches)))
        stop("No way to map colData to ExperimentList")
    else if (!length(matches) && !isEmpty(experiments))
        stop("colData rownames and ExperimentList colnames are empty")
    primary <- rownames(colData)[matches]
    autoMap <- S4Vectors::DataFrame(
        assay=assay, primary=primary, colname=colname)
    missingPrimary <- is.na(autoMap[["primary"]])
    if (nrow(autoMap) && any(missingPrimary)) {
        notFound <- autoMap[missingPrimary, ]
        warning("Data dropped from ExperimentList (element - column):",
            BiocBaseUtils::selectSome(
                paste("\n", notFound[["assay"]], "-", notFound[["colname"]]),
            ), "\nUnable to map to rows of colData", call. = FALSE)
        autoMap <- autoMap[!missingPrimary, ]
    }
    autoMap
}

.sampleMapFromExisting <- function(colData, sampleMap, colnames) {
    colmatches <- match(colnames, sampleMap[["colname"]])
    notFounds <- which(is.na(colmatches))
    if (!length(colmatches) || all(is.na(colmatches))) {
        pmatch <- match(rownames(colData), colnames)
        noPrim <- which(is.na(pmatch))
        pnames <- colnames[na.omit(pmatch)]
        if (!length(pnames) || all(is.na(pnames)))
            stop("No 'colnames' in experiments could be matched:\n  ",
                 paste(BiocBaseUtils::selectSome(colnames), collapse = ", "))
        else if (length(noPrim))
            .warnnomatch("primary", colnames[noPrim])
        DataFrame(primary = pnames, colname = pnames)
    } else {
        if (length(notFounds))
            .warnnomatch("colnames", colmatches[notFounds])
        sampleMap[na.omit(colmatches), c("primary", "colname")]
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @describeIn ExperimentList Get the dimension names for
#' an \code{ExperimentList} using \code{\linkS4class{CharacterList}}
#' @importFrom methods as
#' @export
setMethod("dimnames", "ExperimentList", function(x) {
    list(
        as(lapply(x, function(g) dimnames(g)[[1]]), "CompressedCharacterList"),
        as(lapply(x, function(g) dimnames(g)[[2]]), "CompressedCharacterList")
    )
})

#' @describeIn ExperimentList Get the column names for an \code{ExperimentList}
#'   as a \code{\linkS4class{CharacterList}} slightly more efficiently
#'
#' @importFrom BiocGenerics colnames
#' @inheritParams BiocGenerics::colnames
#'
#' @export
setMethod("colnames", "ExperimentList",
    function(x, do.NULL = TRUE, prefix = "col") {
        as(lapply(x, colnames), "CompressedCharacterList")
    }
)

#' @describeIn ExperimentList Get the row names for an \code{ExperimentList}
#'   as a \code{\linkS4class{CharacterList}} slightly more efficiently
#'
#' @importFrom BiocGenerics rownames
#' @export
setMethod("rownames", "ExperimentList",
    function(x, do.NULL = TRUE, prefix = "row") {
        as(lapply(x, rownames), "CompressedCharacterList")
    }
)

#' @describeIn MultiAssayExperiment Get the dimension names
#' for a \code{MultiAssayExperiment} object
setMethod("dimnames", "MultiAssayExperiment", function(x) {
    dimnames(experiments(x))
})

#' @export
.DollarNames.MultiAssayExperiment <- function(x, pattern = "")
    grep(pattern, names(colData(x)), value = TRUE)

#' @aliases $,MultiAssayExperiment-method
#' @exportMethod $
#' @rdname MultiAssayExperiment-methods
setMethod("$", "MultiAssayExperiment", function(x, name) {
    colData(x)[[name]]
})

#' @importFrom S4Vectors metadata
#' @rdname MultiAssayExperiment-methods
#' @exportMethod metadata
setMethod("metadata", "MultiAssayExperiment", function(x, ...) {
    callNextMethod()
})

#' @importFrom S4Vectors metadata<-
#' @rdname MultiAssayExperiment-methods
#' @exportMethod metadata<-
setReplaceMethod("metadata", c("MultiAssayExperiment", "ANY"),
    function(x, ..., value) {
        callNextMethod()
    }
)

.splitArgs <- function(args) {
    assayArgNames <- c("mcolname", "background", "type",
                       "make.names", "ranges")
    assayArgs <- args[assayArgNames]
    altArgs <- args[!names(args) %in% assayArgNames]
    assayArgs <- Filter(function(x) !is.null(x), assayArgs)
    list(assayArgs, altArgs)
}

.mergeMAE <- function(x, y) {
    if (any(names(x) %in% names(y)))
        stop("Provide unique experiment names")
    expz <- c(experiments(x), experiments(y))
    sampz <- rbind(sampleMap(x), sampleMap(y))
    coldx <- colData(x)
    coldy <- colData(y)
    cdatz <- S4Vectors::merge(x = coldx, y = coldy,
        by = c("row.names", intersect(names(coldx), names(coldy))),
        all = TRUE, sort = FALSE, stringsAsFactors = FALSE)
    rownames(cdatz) <- cdatz[["Row.names"]]
    cdatz <- cdatz[, names(cdatz) != "Row.names", drop = FALSE]
    metaz <- c(metadata(x), metadata(y))
    new("MultiAssayExperiment",
        ExperimentList = expz,
        colData = cdatz,
        sampleMap = sampz,
        metadata = metaz)
}

#' @describeIn MultiAssayExperiment Add a supported data class to the
#' \code{ExperimentList}
#'
#' @param sampleMap \code{c} method: a \code{sampleMap} \code{list} or
#' \code{DataFrame} to guide merge
#' @param mapFrom Either a \code{logical}, \code{character}, or \code{integer}
#' vector indicating the experiment(s) that have an identical colname order as
#' the experiment input(s). If using a character input, the name must match
#' exactly.
#'
#' @examples
#' example("MultiAssayExperiment")
#'
#' ## Add an experiment
#' test1 <- mae[[1L]]
#' colnames(test1) <- rownames(colData(mae))
#'
#' ## Combine current MultiAssayExperiment with additional experiment
#' ## (no sampleMap)
#' c(mae, newExperiment = test1)
#'
#' test2 <- mae[[3L]]
#' c(mae, newExp = test2, mapFrom = 3L)
#'
#' ## Add experiment using experiment name in mapFrom
#' c(mae, RNASeqGeneV2 = test2, mapFrom = "RNASeqGene")
#'
setMethod("c", "MultiAssayExperiment",
    function(x, ..., sampleMap = NULL, mapFrom = NULL) {
    args <- list(...)
    if (!length(args))
        stop("Provide experiments or a 'MultiAssayExperiment' to concatenate")

    if (
        !is.null(sampleMap) && !is.list(sampleMap) &&
        !is(sampleMap, "DataFrame")
    )
        stop(
            "'sampleMap' must be a 'DataFrame', 'data.frame', 'list' or 'NULL'"
        )

    if (identical(length(args), 1L)) {
        input <- args[[1L]]
        if (is(input, "MultiAssayExperiment"))
            return(.mergeMAE(x, input))
        else if (inherits(input, "List") && !is(input, "DataFrame")
                 || is(input, "list"))
            args <- input
    }

    exps <- as(args, "ExperimentList")

    xmap <- sampleMap(x)
    cdata <- colData(x)
    if (!isEmpty(exps)) {
        if (!is.null(mapFrom)) {
            warning("Assuming column order in the data provided ",
                "\n matches the order in 'mapFrom' experiment(s) colnames",
                    call. = FALSE)
            addMaps <- mapToList(sampleMap(x))[mapFrom]
            names(addMaps) <- names(exps)
            sampleMap <- mapply(function(x, y) {
                x[["colname"]] <- colnames(y)
                x
            }, x = addMaps, y = exps)
        } else if (is.null(sampleMap) || isEmpty(sampleMap)) {
            sampleMap <- lapply(colnames(exps), .sampleMapFromExisting,
                colData = cdata, sampleMap = xmap)
        }
        if (!is(sampleMap, "DataFrame"))
            sampleMap <- listToMap(sampleMap)
        newListMap <- rbind(xmap, sampleMap)
        x <- BiocBaseUtils::setSlots(x,
            ExperimentList = c(experiments(x), exps),
            sampleMap = newListMap
        )
    }
    return(x)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Converter - exportClass
###

.metasize <- function(metlist) {
    atmos <- vapply(metlist, is.atomic, logical(1L))
    sum(any(atmos), !atmos)
}

.chr2fxn <- function(fmt) {
    sep <- switch(fmt, csv = ",", "\t")
    cols <- switch(fmt, csv = NA, TRUE)
    qme <- switch(fmt, csv = "double", "escape")
    function(...) utils::write.table(..., sep = sep,
        col.names = cols, qmethod = qme)
}

.sortMetadata <- function(object, objname, dir, fmt, ext, ...) {
    metas <- metadata(object)
    stopifnot(is.list(metas))
    atmos <- vapply(metas, is.atomic, logical(1L))
    metatxt <- metas[atmos]

    fnames <- character(0L)
    if (length(metatxt)) {
        fpath <- file.path(dir, paste0(objname, "_META_0", ext))
        if (!is.function(fmt))
            fmt <- .chr2fxn(fmt)
        fmt(as.data.frame(metatxt), fpath, ...)
        fnames <- fpath
    }
    if (any(!atmos)) {
        nonato <- metas[!atmos]
        nonatos <- seq_along(nonato)
        mpaths <- file.path(dir, paste0(objname, "_META_", nonatos, ext))
        tryCatch({
            invisible(
                lapply(nonatos, function(i) {
                    fmt(as(nonato[[i]], "data.frame"), mpaths[[i]], ...)
                })
            )
        }, error = function(e) conditionMessage(e))
        fnames <- c(fnames, mpaths)
    }
    fnames
}


#' @export
setGeneric(
    "exportClass",
    function(
        object, dir = tempdir(), fmt, ext, match = FALSE, verbose = TRUE, ...
    )
        standardGeneric("exportClass")
)

#' @describeIn MultiAssayExperiment Export data from class to a series
#'     of text files
#'
#' @param dir character(1) A directory for saving exported data (default:
#'     `tempdir()`)
#'
#' @param fmt character(1) or function() Either a format character atomic as
#'     supported by `write.table` either ('csv', or 'tsv') or a function whose
#'     first two arguments are 'object to save' and 'file location'
#'
#' @param ext character(1) A file extension supported by the format argument
#'
#' @param match logical(1) Whether to coerce the current object to a
#'     'MatchedAssayExperiment' object (default: FALSE)
#'
#' @param verbose logical(1) Whether to print additional information (default
#'     TRUE)
#'
#' @aliases exportClass
#' @exportMethod exportClass
setMethod("exportClass", "MultiAssayExperiment",
    function(object, dir = tempdir(), fmt, ext, match = FALSE,
            verbose = TRUE, ...) {
        if (missing(dir) || !dir.exists(dir))
            stop("Specify a valid folder location for saving data files")
        objname <- as.character(substitute(object))

        if (isTRUE(match))
            object <- as(object, "MatchedAssayExperiment")

        nfiles <- sum(length(object), !isEmpty(colData(object)),
            !isEmpty(sampleMap(object)), .metasize(metadata(object)))
        if (missing(verbose) || !isFALSE(verbose))
            message("Writing about ", nfiles, " files to disk...")

        if (is.function(fmt) && missing(ext))
            stop("Provide a valid file extention, see 'ext' argument")

        if (!is.function(fmt) && !is.character(fmt))
            stop("Invalid format type: must be 'character' or 'function'")

        if (is.character(fmt)) {
            if (missing(ext))
                ext <- paste0(".", fmt)
            fmt <- .chr2fxn(fmt)
        }

        coldatname <- paste0(objname, "_", "colData", ext)
        sampmapname <- paste0(objname, "_", "sampleMap", ext)

        exfnames <- file.path(dir,
            c(paste0(objname, "_", names(experiments(object)), ext),
            coldatname, sampmapname))
        alists <- lapply(assays(object), as.data.frame)
        lists <- c(alists, list(coldat = as.data.frame(colData(object)),
            sampmap = as.data.frame(sampleMap(object))))

        metafs <- .sortMetadata(object, objname, dir, fmt, ext)

        invisible(Map(function(fname, lobject) {
            fmt(lobject, fname, ...)
        }, exfnames, lists))

        c(metafs, exfnames)
    }
)

