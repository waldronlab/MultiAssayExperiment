#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' IRanges
NULL

## Helper function for validity checks
.uniqueSortIdentical <- function(charvec1, charvec2) {
    listInput <- list(charvec1, charvec2)
    listInput <- lapply(listInput, function(x) sort(unique(x)))
    return(identical(listInput[[1]], listInput[[2]]))
}

.allIn <- function(charvec1, charvec2) {
    return(all(charvec2 %in% charvec1))
}

### ==============================================
### MultiAssayExperiment class
### ----------------------------------------------

#' MultiAssayExperiment - An integrative multi-assay class for experiment data
#'
#' @description
#' The `MultiAssayExperiment` class can be used to manage results of
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of
#' \code{\linkS4class{SummarizedExperiment}},
#' \code{\linkS4class{ExpressionSet}}, `matrix`,
#' \code{\link[RaggedExperiment:RaggedExperiment-class]{RaggedExperiment}}
#' (inherits from \code{\linkS4class{GRangesList}}), and `RangedVcfStack`.
#' Create new `MultiAssayExperiment` instances with the homonymous
#' constructor, minimally with the argument \code{\link{ExperimentList}},
#' potentially also with the arguments `colData` (see section below) and
#' \code{\link{sampleMap}}.
#'
#' @details
#' The dots (\code{\ldots}) argument allows the user to specify additional
#' arguments in several instances.
#' \itemize{
#'     \item subsetting \strong{[}: additional arguments sent to
#'         \link[GenomicRanges:findOverlaps-methods]{findOverlaps}.
#'     \item mergeReplicates: used to specify arguments for the \code{simplify}
#'         functional argument
#'     \item assay: may contain withDimnames, which is forwarded to assays
#'     \item combining \strong{c}: compatible \code{MultiAssayExperiment} classes
#'         passed on to the \code{\linkS4class{ExperimentList}} constructor,
#'         can be a \code{list}, \code{\linkS4class{List}}, or a series of
#'         named arguments. See the examples below.
#' }
#'
#' @section colData:
#' The `colData` slot is a collection of primary specimen data valid
#' across all experiments. This slot is strictly of class
#' \code{\linkS4class{DataFrame}} but arguments for the constructor function
#' allow arguments to be of class `data.frame` and subsequently coerced.
#'
#' @section ExperimentList:
#' The \code{\link{ExperimentList}} slot is designed to contain results from
#' each experiment/assay. It contains a \linkS4class{SimpleList}.
#'
#' @section sampleMap:
#' The \code{\link{sampleMap}} contains a `DataFrame` of translatable
#' identifiers of samples and participants or biological units. The standard
#' column names of the `sampleMap` are "assay", "primary", and "colname".
#' Note that the "assay" column is a factor corresponding to the names of each
#' experiment in the `ExperimentList`. In the case where these names do
#' not match between the `sampleMap` and the experiments, the documented
#' experiments in the `sampleMap` take precedence and experiments are
#' dropped by the harmonization procedure. The constructor function will
#' generate a `sampleMap` in the case where it is not provided and this
#' method may be a 'safer' alternative for creating the `MultiAssayExperiment`
#' (so long as the rownames are identical in the `colData`, if provided).
#' An empty `sampleMap` may produce empty experiments if the levels of the
#' "assay" factor in the `sampleMap` do not match the names in the
#' `ExperimentList`.
#'
#' @slot ExperimentList A \code{\link{ExperimentList}} class object for
#' each assay dataset
#'
#' @slot colData A `DataFrame` of all clinical/specimen data available
#' across experiments
#'
#' @slot sampleMap A `DataFrame` of translatable identifiers
#' of samples and participants
#'
#' @slot metadata Additional data describing the
#' `MultiAssayExperiment` object
#'
#' @slot drops A metadata `list` of dropped information
#'
#' @param object,x A `MultiAssayExperiment` object
#'
#' @param ... Additional arguments for supporting functions. See details.
#'
#' @return A `MultiAssayExperiment` object
#'
#' @md
#'
#' @examples
#' example("MultiAssayExperiment")
#'
#' ## Subsetting
#' # Rows (i) Rows/Features in each experiment
#' mae[1, , ]
#' mae[c(TRUE, FALSE), , ]
#'
#' # Columns (j) Rows in colData
#' mae[, rownames(colData(mae))[3:2],  ]
#'
#' # Assays (k)
#' mae[, , "Affy"]
#'
#' ## Complete cases (returns logical vector)
#' completes <- complete.cases(mae)
#' compMAE <- mae[, completes, ]
#' compMAE
#' colData(compMAE)
#'
#' @exportClass MultiAssayExperiment
#'
#' @seealso \link{MultiAssayExperiment-methods} for slot modifying methods,
#'     \href{https://github.com/waldronlab/MultiAssayExperiment/wiki/MultiAssayExperiment-API}{MultiAssayExperiment API}
#'
#' @include ExperimentList-class.R
setClass(
    "MultiAssayExperiment",
    contains = "Annotated",
    slots = list(
        ExperimentList = "ExperimentList",
        colData = "DataFrame",
        sampleMap = "DataFrame",
        drops = "list"
    ),
    prototype = prototype(
        colData = new("DFrame"),
        sampleMap = new("DFrame")
    )
)

### ==============================================
### MultiAssayExperiment constructor
### ----------------------------------------------

.harmonize <- function(experiments, colData, sampleMap) {
    harmony <- character()
    ## sampleMap assays agree with experiment names
    assay <- intersect(levels(sampleMap[["assay"]]), names(experiments))
    keep_sampleMap_assay <- sampleMap[["assay"]] %in% assay
    if (!all(keep_sampleMap_assay)) {
        sampleMap <- sampleMap[keep_sampleMap_assay, , drop=FALSE]
        sampleMap[["assay"]] <- factor(sampleMap[["assay"]], levels=assay)
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_assay),
                  "sampleMap rows not in names(experiments)"))
    }

    ## sampleMap colname agrees with experiment colnames
    grp <- droplevels(sampleMap[["assay"]])
    colnm <- split(sampleMap[["colname"]], grp)
    keep <- Map(intersect, colnm, colnames(experiments)[names(colnm)])
    keep_sampleMap_colname <- logical(nrow(sampleMap))
    split(keep_sampleMap_colname, grp) <- Map("%in%", colnm, keep)
    if (!all(keep_sampleMap_colname)) {
        sampleMap <- sampleMap[keep_sampleMap_colname, , drop=FALSE]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_colname),
                  "sampleMap rows with 'colname'",
                  "not in colnames of experiments"))
    }

    ## sampleMap primary agrees with primary
    primary <- intersect(rownames(colData), sampleMap[["primary"]])
    keep_sampleMap_primary <- sampleMap[["primary"]] %in% primary
    if (!all(keep_sampleMap_primary)) {
        sampleMap <- sampleMap[keep_sampleMap_primary, , drop=FALSE]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_sampleMap_primary),
                  "sampleMap rows with 'primary' not in colData"))
    }

    ## update objects
    assay <- intersect(names(experiments), levels(sampleMap[["assay"]]))
    experiments_columns <- split(sampleMap[["colname"]], sampleMap[["assay"]])
    primary <- intersect(rownames(colData), sampleMap[["primary"]])
    keep_colData <- rownames(colData) %in% primary
    if (!all(keep_colData)) {
        colData <- colData[keep_colData, , drop = FALSE]
        harmony <- c(
            harmony,
            paste("removing", sum(!keep_colData),
                  "colData rownames not in sampleMap 'primary'"))
    }

    experiments <- mendoapply(function(x, idx) {
        colmatch <- colnames(x) %in% idx
        if (!all(colmatch))
            x <- x[, colmatch, drop=FALSE]
        x
    }, experiments[assay], experiments_columns[assay])

    if (length(harmony))
        message("harmonizing input:\n  ", paste(harmony, collapse="\n  "))
    list(experiments=experiments, sampleMap=sampleMap, colData=colData)
}

.checkFixSampleMap <- function(samplemap) {
    samplemap <- as(samplemap, "DataFrame")
    smapnames <- c("assay", "primary", "colname")
    if (!all(smapnames %in% colnames(samplemap)))
        stop("'sampleMap' does not have required columns")

    if (!identical(colnames(samplemap), smapnames))
        samplemap <- samplemap[, smapnames]

    .smapColumnCoerce(samplemap)
}

.smapColumnCoerce <- function(samplemap) {
    isfuns <- list(is.factor, is.character, is.character)
    colsOK <- mapply(
        FUN = function(fun, X) { fun(X) },  fun = isfuns, X = samplemap
    )
    if (all(colsOK))
        return(samplemap)
    asfuns <- list(
        as.factor = as.factor,
        as.character = as.character,
        as.character = as.character
    )
    samplemap[] <- Map(
        function(cName, isFun, coerceFun, funname) {
            smapCol <- samplemap[[cName]]
            if (!isFun(smapCol))
                warning(
                    "sampleMap[['", cName, "']] coerced with ", funname, "()",
                    call. = FALSE
                )
            samplemap[[cName]] <- coerceFun(samplemap[[cName]])
        },
        cName = names(samplemap),
        isFun = isfuns,
        coerceFun = asfuns,
        funname = names(asfuns)
    )

    samplemap
}

#' Construct an integrative representation of multi-omic data with
#' \code{MultiAssayExperiment}
#'
#' The constructor function for the \link{MultiAssayExperiment-class} combines
#' multiple data elements from the different hierarchies of data
#' (study, experiments, and samples). It can create instances where neither
#' a \code{sampleMap} or a \code{colData} set is provided. Please see the
#' MultiAssayExperiment API documentation for more information.
#'
#' @section colData:
#' The `colData` input can be either `DataFrame` or `data.frame` with
#' subsequent coercion to DataFrame. The rownames in the `colData` must match
#' the colnames in the experiments if no sampleMap is provided.
#'
#' @section experiments:
#' The `experiments` input can be of class \linkS4class{SimpleList} or `list`.
#' This input becomes the \code{\link{ExperimentList}}. Each element of the
#' input `list` or `List` must be named, rectangular with two dimensions, and
#' have `dimnames`.
#'
#' @section sampleMap:
#' The \code{\link{sampleMap}} can either be input as `DataFrame` or
#' `data.frame` with eventual coercion to `DataFrame`. The `sampleMap` relates
#' biological units and biological measurements within each assay. Each row in
#' the `sampleMap` is a single such link. The standard column names of the
#' `sampleMap` are "assay", "primary", and "colname".  Note that the "assay"
#' column is a factor corresponding to the names of each experiment in the
#' `ExperimentList`. In the case where these names do not match between the
#' `sampleMap` and the experiments, the documented experiments in the
#' `sampleMap` take precedence and experiments are dropped by the harmonization
#' procedure. The constructor function will generate a `sampleMap` in the case
#' where it is not provided and this method may be a 'safer' alternative for
#' creating the `MultiAssayExperiment` (so long as the rownames are identical
#' in the `colData`, if provided).  An empty `sampleMap` may produce empty
#' experiments if the levels of the "assay" factor in the `sampleMap` do not
#' match the names in the `ExperimentList`.
#'
#' @param experiments A \code{list} or \link{ExperimentList} of all
#' combined experiments
#'
#' @param colData A \code{\linkS4class{DataFrame}} or \code{data.frame} of
#' characteristics for all biological units
#'
#' @param sampleMap A \code{DataFrame} or \code{data.frame} of assay names,
#' sample identifiers, and colname samples
#'
#' @param metadata An optional argument of "ANY" class (usually list) for
#' content describing the experiments
#'
#' @param drops A \code{list} of unmatched information
#' (included after subsetting)
#'
#' @return A \code{MultiAssayExperiment} object that can store
#' experiment and phenotype data
#'
#' @example inst/scripts/MultiAssayExperiment-Ex.R
#'
#' @export MultiAssayExperiment
#' @seealso \link{MultiAssayExperiment-class}
MultiAssayExperiment <-
    function(
        experiments = ExperimentList(),
        colData = S4Vectors::DataFrame(),
        sampleMap = S4Vectors::DataFrame(
            assay = factor(),
            primary = character(),
            colname = character()
        ),
        metadata = list(),
        drops = list()
    )
{
    experiments <- as(experiments, "ExperimentList")

    allsamps <- unique(unlist(unname(colnames(experiments))))
    if (missing(sampleMap)) {
        if (missing(colData))
            colData <- S4Vectors::DataFrame(row.names = allsamps)
        sampleMap <- .sampleMapFromData(colData, experiments)
    } else {
        if (missing(colData))
            colData <- S4Vectors::DataFrame(
                row.names = unique(sampleMap[["primary"]])
            )
    }

    sampleMap <- .checkFixSampleMap(sampleMap)

    colData <- as(colData, "DataFrame")

    bliss <- .harmonize(experiments, colData, sampleMap)

    new("MultiAssayExperiment",
        ExperimentList = bliss[["experiments"]],
        colData = bliss[["colData"]],
        sampleMap = bliss[["sampleMap"]],
        metadata = metadata)
}

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

## EXPERIMENTLIST
## 1.i. ExperimentList length must be the same as the unique length of the
## sampleMap "assay" column.
.checkExperimentList <- function(object) {
    errors <- character()
    assays <- levels(sampleMap(object)[["assay"]])
    if (!is.null(assays) &&
        !identical(length(experiments(object)), length(assays))) {
        msg <- paste0("ExperimentList must be the same length as",
                      " the sampleMap assay column")
        errors <- c(errors, msg)
    }

## 1.ii. Element names of the ExperimentList should be found in the
## sampleMap "assay" column.
    if (!is.null(assays) && !all(names(experiments(object)) %in% assays)) {
        msg <- paste0("All ExperimentList names were not found in",
                      " the sampleMap assay column")
        errors <- c(errors, msg)
    }

## 1.iii. ExperimentList cannot have any non-empty elements when the sampleMap
## is empty.
    if (isEmpty(sampleMap(object)) && !isEmpty(experiments(object))) {
        msg <- paste0(
            "All non-empty elements in the ExperimentList must have",
            " names in the sampleMap assay column"
        )
        errors <- c(errors, msg)
    }

    if (!length(errors)) NULL else errors
}

.id_sort_uniq <- function(x, y) {
    identical(
        sort(unique(x)), sort(unique(y))
    )
}

## 1.iii. For each ExperimentList element, colnames must be found in the
## "assay" column of the sampleMap
.checkSampleNames <- function(object) {
    sampMap <- sampleMap(object)
    assayCols <- mapToList(sampMap[, c("assay", "colname")])
    colNams <- Filter(function(x) !isEmpty(x), colnames(object))
    msg <- NULL
    if (length(colNams)) {
        logicResult <- mapply(function(x, y) {
            .id_sort_uniq(x = x, y = y)
        }, x = colNams, y = assayCols[names(colNams)])
        if (!all(logicResult))
            msg <- "not all ExperimentList samples are found in the sampleMap"
    }
    msg
}

## COLDATA
## 2.i. See setClass above where colData = "DataFrame"

## SAMPLEMAP
## 3.i. all values in the sampleMap "primary" column must be found in the
## rownames of colData
## 3.i.a sampleMap assay column must be a factor
.checkSampleMapNamesClass <- function(object) {
    errors <- character()
    if (!(.allIn(
        rownames(colData(object)),
        sampleMap(object)[["primary"]]
    ))) {
        msg <- "All 'sampleMap[[primary]]' must be in 'rownames(colData)'"
        errors <- c(errors, msg)
    }
    if (!is.factor(sampleMap(object)[["assay"]])) {
        msg <- "'sampleMap' assay column not a factor"
        errors <- c(errors, msg)
    }
    if (!length(errors)) NULL else errors
}

## 3.ii. Within rows of "sampleMap" corresponding to a single value in the
## "assay" column, there can be no duplicated values in the "colname" column
.uniqueNamesInAssays <- function(object) {
    lcheckdups <- colnames(object)
    logchecks <- any(vapply(lcheckdups, FUN = function(x) {
        as.logical(anyDuplicated(x))
    }, FUN.VALUE = logical(1L)))
    if (!logchecks)
        NULL
    else
        "All colname identifiers in assays must be unique"
}

.validMultiAssayExperiment <- function(object) {
    if (length(experiments(object)) != 0L) {
        c(.checkExperimentList(object),
          .checkSampleMapNamesClass(object),
          .uniqueNamesInAssays(object),
          .checkSampleNames(object)
        )
    }
}

S4Vectors::setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)

.hasOldAPI <- function(object) {
    isTRUE(.hasSlot(object, "Elist")) || isTRUE(.hasSlot(object, "pData"))
}

#' @exportMethod show
#' @describeIn MultiAssayExperiment Show method for a
#' \code{MultiAssayExperiment}
setMethod("show", "MultiAssayExperiment", function(object) {
    if (.hasOldAPI(object)) {
        stop("MultiAssayExperiment is outdated, please run updateObject()")
    }
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    if (!length(o_names)) {
        o_names <- "none"
    }
    c_elist <- class(experiments(object))
    c_mp <- class(colData(object))
    c_sm <- class(sampleMap(object))
    cat(sprintf("A %s", o_class),
        "object of", o_len, "listed\n",
        ifelse(o_len == 1L, "experiment", "experiments"),
        "with",
        ifelse(identical(o_names, "none"), "no user-defined names",
               ifelse(length(o_names) == 1L, "a user-defined name",
                      "user-defined names")),
        ifelse(length(o_len) == 0L, "or", "and"),
        ifelse(length(o_len) == 0L, "classes.",
               ifelse(o_len == 1L,
                      "respective class.\n", "respective classes.\n")),
        "Containing an ")
    show(experiments(object))
    cat("Functionality:\n experiments() - obtain the ",
        sprintf("%s", c_elist), " instance",
        "\n colData() - the primary/phenotype DataFrame",
        "\n sampleMap() - the sample coordination DataFrame",
        "\n `$`, `[`, `[[` - extract colData columns, subset, or experiment",
        "\n *Format() - convert ", "into a long or wide DataFrame",
        "\n assays() - convert ", sprintf("%s", c_elist),
        " to a SimpleList of matrices",
        "\n exportClass() - save data to flat files\n",
    sep = "")
})

#' @name MultiAssayExperiment-methods
#' @title Accessing and modifying information in MultiAssayExperiment
#'
#' @description A set of accessor and setter generic functions to extract
#' either the \code{sampleMap}, the \code{\link{ExperimentList}},
#' \code{colData}, or \code{metadata} slots of a
#' \code{\link{MultiAssayExperiment}} object
#'
#' @section Accessors:
#' Eponymous names for accessing \code{MultiAssayExperiment} slots with the
#' exception of the \link{ExperimentList} accessor named \code{experiments}.
#' \itemize{
#'    \item colData: Access the \code{colData} slot
#'    \item sampleMap: Access the \code{sampleMap} slot
#'    \item experiments: Access the \link{ExperimentList} slot
#'    \item `[[`: Access the \link{ExperimentList} slot
#'    \item `$`: Access a column in \code{colData}
#'    \item `drops`: Get a vector of dropped \link{ExperimentList} names
#' }
#'
#' @section Setters:
#' Setter method values (i.e., '\code{function(x) <- value}'):
#' \itemize{
#'     \item experiments<-: An \code{\link{ExperimentList}} object
#'     containing experiment data of supported classes
#'     \item sampleMap<-: A \code{\link{DataFrame}} object relating
#'     samples to biological units and assays
#'     \item colData<-: A \code{\link{DataFrame}} object describing the
#'     biological units
#'     \item metadata<-: A \code{list} object of metadata
#'     \item `[[<-`: Equivalent to the \code{experiments<-} setter method for
#'     convenience
#'     \item `$<-`: A vector to replace the indicated column in \code{colData}
#'     \item `drops<-`: Trace \link{ExperimentList} names that have been
#'     removed
#' }
#'
#' @param object,x A \code{MultiAssayExperiment} object
#'
#' @param name A column in \code{colData}
#'
#' @param value See details.
#'
#' @param ... Argument not in use
#'
#' @return Accessors: Either a \code{sampleMap}, \code{ExperimentList}, or
#' \code{DataFrame} object
#' @return Setters: A \code{MultiAssayExperiment} object
#'
#' @example inst/scripts/MultiAssayExperiment-methods-Ex.R
#'
#' @aliases experiments sampleMap experiments<- sampleMap<- drops drops<-
NULL

### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

#' @export
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))

#' @exportMethod sampleMap
#' @rdname MultiAssayExperiment-methods
setMethod("sampleMap", "MultiAssayExperiment", function(x)
    getElement(x, "sampleMap"))

#' @export
setGeneric("experiments", function(x) standardGeneric("experiments"))

#' @exportMethod experiments
#' @rdname MultiAssayExperiment-methods
setMethod("experiments", "MultiAssayExperiment", function(x)
    getElement(x, "ExperimentList"))


#' @exportMethod colData
#' @rdname MultiAssayExperiment-methods
setMethod("colData", "MultiAssayExperiment", function(x, ...) {
    getElement(x, "colData")
})

#' @export
setGeneric("drops", function(x) standardGeneric("drops"))

#' @exportMethod drops
#' @rdname MultiAssayExperiment-methods
setMethod("drops", "MultiAssayExperiment", function(x) {
    getElement(x, "drops")
})

### - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @exportMethod length
#' @describeIn MultiAssayExperiment Get the length of ExperimentList
setMethod("length", "MultiAssayExperiment", function(x)
    length(experiments(x))
)

#' @exportMethod names
#' @describeIn MultiAssayExperiment Get the names of the ExperimentList
setMethod("names", "MultiAssayExperiment", function(x)
    names(experiments(x))
)

### - - - - - - - - - - - - - - - - - - - - - - - -
### Replacers
###

#' @export
setGeneric("sampleMap<-", function(object, value)
    standardGeneric("sampleMap<-"))

#' @exportMethod sampleMap<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "DataFrame"),
    function(object, value) {
        if (isEmpty(value))
            value <- DataFrame(assay = factor(), primary = character(),
                colname = character())
        rebliss <- .harmonize(experiments(object),
            colData(object),
            value)

        BiocBaseUtils::setSlots(object,
            ExperimentList = rebliss[["experiments"]],
            colData = rebliss[["colData"]],
            sampleMap = rebliss[["sampleMap"]],
            check = FALSE
        )
    }
)

#' @rdname MultiAssayExperiment-methods
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "ANY"),
    function(object, value) {
        stopifnot(is.data.frame(value))
        value <- as(value, "DataFrame")
        `sampleMap<-`(object, value)
    }
)

#' @export
setGeneric("experiments<-", function(object, value)
    standardGeneric("experiments<-"))

#' @export
#' @rdname MultiAssayExperiment-methods
setGeneric("drops<-", function(x, ..., value) standardGeneric("drops<-"))

#' @exportMethod experiments<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("experiments", c("MultiAssayExperiment", "ExperimentList"),
    function(object, value) {
        if (!any(names(object) %in% names(value)) && !isEmpty(object)) {
            if (isEmpty(drops(object)))
                warning("'experiments' dropped; see 'drops()'", call. = FALSE)
            drops(object) <-
                list(experiments = setdiff(names(object), names(value)))
        }
        o_cnames <- colnames(object)
        v_cnames <- colnames(value)
        if (identical(o_cnames, v_cnames)) {
            BiocBaseUtils::setSlots(
                object = object,
                ExperimentList = value,
                check = FALSE
            )
        } else {

            samplemap <- sampleMap(object)

            if (all(names(o_cnames) %in% names(v_cnames))) {
                levels <- names(v_cnames)
                ordernames <- names(Filter(length, v_cnames))
                samplemap <- listToMap(mapToList(samplemap)[ordernames])
                samplemap[["assay"]] <-
                    factor(samplemap[["assay"]], levels = levels)
            }

            rebliss <- .harmonize(value, colData(object), samplemap)

            BiocBaseUtils::setSlots(
                object = object,
                ExperimentList = rebliss[["experiments"]],
                colData = rebliss[["colData"]],
                sampleMap = rebliss[["sampleMap"]],
                check = FALSE
            )
        }
    }
)

#' @exportMethod experiments<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("experiments", c("MultiAssayExperiment", "List"),
    function(object, value) {
        value <- as(value, "ExperimentList")
        experiments(object) <- value
        object
    }
)


#' @exportMethod colData<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("colData", c("MultiAssayExperiment", "DataFrame"),
    function(x, value) {
        x_rnames <- rownames(colData(x))
        v_rnames <- rownames(value)

        if (!any(v_rnames %in% x_rnames) && !isEmpty(value))
            stop("'rownames(value)' have no match in 'rownames(colData)';\n  ",
                "See '?renamePrimary' for renaming primary units")

        if (identical(x_rnames, v_rnames))
            BiocBaseUtils::setSlots(
                object = x,
                colData = value,
                check = FALSE
            )
        else {
            rebliss <- .harmonize(experiments(x), value, sampleMap(x))

            BiocBaseUtils::setSlots(
                object = x,
                ExperimentList = rebliss[["experiments"]],
                colData = rebliss[["colData"]],
                sampleMap = rebliss[["sampleMap"]],
                check = FALSE
            )
        }
    }
)

#' @rdname MultiAssayExperiment-methods
setReplaceMethod("colData", c("MultiAssayExperiment", "ANY"),
    function(x, value) {
        if (!is.data.frame(value))
            stop("'colData' can be either 'data.frame' or 'DataFrame'")
        value <- as(value, "DataFrame")
        `colData<-`(x, value)
    }
)

#' @exportMethod drops<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("drops", c("MultiAssayExperiment", "ANY"),
    function(x, ..., value) {
        anydrops <- getElement(x, "drops")[["experiments"]]
        if (length(anydrops))
            value[["experiments"]] <- union(anydrops, value[["experiments"]])
        slot(x, "drops") <- value
        return(x)
})

#' @exportMethod $<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("$", "MultiAssayExperiment", function(x, name, value) {
    colData(x)[[name]] <- value
    return(x)
})

#' @exportMethod names<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("names", c("MultiAssayExperiment", "ANY"),
    function(x, value)
{
    if (!is.character(value))
        stop("'value' must be a character vector",
                "in names(x) <- value")
    if (length(value) != length(x))
        stop("experiment names and experiments not equal in length")

    explist <- experiments(x)
    oldNames <- names(explist)
    names(explist) <- value
    sampmap <- sampleMap(x)
    map <- .setNames(value, oldNames)
    sampmap[, "assay"] <-
        factor(unname(map[sampmap[["assay"]]]), levels = value)

    BiocBaseUtils::setSlots(x,
        ExperimentList = explist,
        sampleMap = sampmap,
        check = FALSE)
})

#' @exportMethod colnames<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("colnames", c("MultiAssayExperiment", "List"),
    function(x, value)
{
    if (!is(value, "CharacterList"))
        stop("'value' must be a 'CharacterList' in 'colnames(x) <- value'")
    if (length(value) != length(x))
        stop("'colnames(x)' and 'value' not equal in length")

    cnames <- colnames(x)
    if (!identical(lengths(value), lengths(cnames)))
        stop("'value' names and lengths should all be identical to 'names(x)'")

    samplemap <- sampleMap(x)
    splitmap <- mapToList(samplemap)
    splitmap <- Map(function(x, y) {
        x[["colname"]] <- y
        x
    }, x = splitmap, y = value)
    exps <- experiments(x)
    exps <- S4Vectors::mendoapply(function(x, y) {
        colnames(x) <- y
        x
    }, x = exps, y = value)

    BiocBaseUtils::setSlots(
        object = x,
        ExperimentList = exps,
        sampleMap = listToMap(splitmap)
    )
})

#' @exportMethod colnames<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("colnames", c("MultiAssayExperiment", "list"),
    function(x, value)
{
    value <- as(value, "CharacterList")
    colnames(x) <- value
    x
})


#' @exportMethod updateObject
#'
#' @param verbose logical (default FALSE) whether to print extra messages
#'
#' @describeIn MultiAssayExperiment Update old serialized MultiAssayExperiment
#' objects to new API
setMethod("updateObject", "MultiAssayExperiment",
    function(object, ..., verbose = FALSE) {
        if (verbose)
            message("updateObject(object = 'MultiAssayExperiment')")

        oldEL <- try(object@ExperimentList, silent = TRUE)
        if (is(oldEL, "try-error")) {
            explist <- ExperimentList(object@Elist@listData)
            samplemap <- .checkFixSampleMap(object@sampleMap)
        } else {
            explist <- experiments(object)
            samplemap <- sampleMap(object)
        }

        oldCD <- try(object@colData, silent = TRUE)
        if (is(oldCD, "try-error"))
            coldata <- object@pData
        else
            coldata <- colData(object)

        explist <- updateObject(explist, ..., verbose = verbose)
        coldata <- updateObject(coldata, ..., verbose = verbose)
        samplemap <- updateObject(samplemap, ..., verbose = verbose)

        BiocBaseUtils::setSlots(
            object,
            ExperimentList = explist,
            colData = coldata,
            sampleMap = samplemap,
            check=FALSE
        )
    }
)

.mergeColData <- function(inlist) {
    CDbyEXP <- lapply(names(inlist),
        function(i, x) {
            tryCatch({
                S4Vectors::DataFrame(colData(x[[i]]), experiment_name = i)
            }, error = function(e) {
                S4Vectors::DataFrame(row.names = colnames(x[[i]]))
            } )
        }, x = inlist
    )
    colDatas <- Filter(function(y) { !isEmpty(y) }, CDbyEXP)
    if (length(colDatas)) {
        rnames <- unlist(lapply(colDatas, rownames))
        res <- Reduce(function(x, y) {
            S4Vectors::merge(
                x, y, by = intersect(names(x), names(y)),
                all = TRUE, sort = FALSE
            )
        }, colDatas)
        rownames(res) <- rnames
    } else {
        res <- S4Vectors::DataFrame(
            row.names = unlist(lapply(CDbyEXP, rownames))
        )
    }
    res
}

#' @rdname MultiAssayExperiment-class
#'
#' @name coerce-MultiAssayExperiment
#'
#' @aliases coerce,list,MultiAssayExperiment-method
#'     coerce,List,MultiAssayExperiment-method
#'
#' @section
#' coercion:
#'   Convert a `list` or S4 `List` to a MultiAssayExperiment object using the
#'   \link[methods]{as} function.
#'
#' In the following example, `x` is either a `list` or \linkS4class{List}:
#'
#'   `as(x, "MultiAssayExperiment")`
#'
#'   Convert a `MultiAssayExperiment` to `MAF` class object using the
#'   \link[methods]{as} function.
#'
#' In the following example, `x` is a \linkS4class{MultiAssayExperiment}:
#'
#'   `MultiAssayExperimentToMAF(x)`
#'
#' @md
#'
#' @exportMethod coerce

setAs("list", "MultiAssayExperiment", function(from) {
        newfrom <- as(from, "List")
        as(newfrom, "MultiAssayExperiment")
    }
)

setAs("List", "MultiAssayExperiment", function(from) {
        metaf <- metadata(from)
        explist <- as(from, "ExperimentList")
        colData <- .mergeColData(from)
        MultiAssayExperiment(
            experiments = explist, colData = colData, metadata = metaf
        )
    }
)

