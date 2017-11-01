#' @import BiocGenerics SummarizedExperiment S4Vectors GenomicRanges methods
#' @importFrom Biobase pData
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom Biobase pData<-
#' @importFrom S4Vectors metadata<-
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

#' An integrative multi-assay class for experiment data
#'
#' @description
#' The \code{MultiAssayExperiment} class can be used to manage results of
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of
#' \code{\linkS4class{SummarizedExperiment}},
#' \code{\linkS4class{ExpressionSet}},
#' \code{matrix}, \code{\link[RaggedExperiment]{RaggedExperiment}}
#' (inherits from \code{\linkS4class{GRangesList}}), and \code{RangedVcfStack}.
#' Create new \code{MultiAssayExperiment} instances with the homonymous
#' constructor, minimally with the argument \code{\link{ExperimentList}},
#' potentially also with the arguments \code{colData} (see section below) and
#' \code{\link{sampleMap}}.
#'
#' @details
#' The dots (\code{\ldots}) argument allows the user to specify additional
#' arguments in several instances. When subsetting (\strong{[}) a
#' \code{MultiAssayExperiment}, the dots allow for additional
#' arguments to be sent to \link{findOverlaps}. When using the
#' \code{mergeReplicates} method, the dots are used to specify arguments for
#' the supplied \code{simplify} argument and function. When using the
#' \strong{assay} method. When using \strong{c} method
#' to add experiments to a \code{MultiAssayExperiment}, the dots allow extra
#' data classes compatible with the MultiAssayExperiment API. See: \link{API}
#'
#' @section colData:
#' The \code{colData} slot is a collection of primary specimen data valid
#' across all experiments. This slot is strictly of class
#' \code{\linkS4class{DataFrame}} but arguments for the constructor function
#' allow arguments to be of class \code{data.frame} and subsequently coerced.
#'
#' @section ExperimentList:
#' The \code{\link{ExperimentList}} slot is designed to contain results from
#' each experiment/assay. It contains a \link[S4Vectors]{SimpleList}.
#'
#' @section sampleMap:
#' The \code{\link{sampleMap}} contains a \code{DataFrame} of translatable
#' identifiers of samples and participants or biological units. Standard column
#' names of the sampleMap are "assay", "primary", and "colname".
#'
#' @slot ExperimentList A \code{\link{ExperimentList}} class object for
#' each assay dataset
#' @slot colData A \code{DataFrame} of all clinical/specimen data available
#' across experiments
#' @slot sampleMap A \code{DataFrame} of translatable identifiers
#' of samples and participants
#' @slot metadata Additional data describing the
#' \code{MultiAssayExperiment} object
#' @slot drops A metadata \code{list} of dropped information
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param ... Additional arguments for supporting functions. See details.
#'
#' @return A \code{MultiAssayExperiment} object
#'
#' @examples
#' example("MultiAssayExperiment")
#'
#' ## Subsetting
#' # Rows (i) Rows/Features in each experiment
#' myMultiAssayExperiment[1, , ]
#' myMultiAssayExperiment[c(TRUE, FALSE), , ]
#'
#' # Columns (j) Rows in colData
#' myMultiAssayExperiment[, rownames(colData(myMultiAssayExperiment))[3:2],  ]
#'
#' # Assays (k)
#' myMultiAssayExperiment[, , "Affy"]
#'
#' ## Complete cases (returns logical vector)
#' completes <- complete.cases(myMultiAssayExperiment)
#' compMAE <- myMultiAssayExperiment[, completes, ]
#' compMAE
#' colData(compMAE)
#'
#' @exportClass MultiAssayExperiment
#' @seealso \link{MultiAssayExperiment-methods} for slot modifying methods
#' @include ExperimentList-class.R
setClass("MultiAssayExperiment",
         slots = list(
           ExperimentList = "ExperimentList",
           colData = "DataFrame",
           sampleMap = "DataFrame",
           metadata = "ANY",
           drops = "list"
         )
)

### - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

## ELIST
## 1.i. ExperimentList length must be the same as the unique length of the
## sampleMap "assay" column.
.checkExperimentList <- function(object) {
    errors <- character()
    assays <- levels(sampleMap(object)[["assay"]])
    if (length(experiments(object)) != length(assays)) {
        msg <- paste0("ExperimentList must be the same length as",
                      " the sampleMap assay column")
        errors <- c(errors, msg)
    }

## 1.ii. Element names of the ExperimentList should be found in the
## sampleMap "assay" column.
    if (!all(names(experiments(object)) %in% assays)) {
        msg <- paste0("All ExperimentList names were not found in",
                      " the sampleMap assay column")
        errors <- c(errors, msg)
    }
    if (length(errors)) NULL else errors
}

## 1.iii. For each ExperimentList element, colnames must be found in the
## "assay" column of the sampleMap
.checkSampleNames <- function(object) {
    sampMap <- sampleMap(object)
    assayCols <- mapToList(sampMap[, c("assay", "colname")])
    colNams <- colnames(object)
    logicResult <- mapply(function(columnNames, assayColumns) {
        identical(sort(columnNames), sort(assayColumns))
    }, columnNames = colNams,
    assayColumns = assayCols)
    if (!all(logicResult)) {
        "not all samples in the ExperimentList are found in the sampleMap"
    }
    NULL
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
        msg <- "All samples in the 'sampleMap' must be in the 'colData'"
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
    SampMap <- sampleMap(object)
    lcheckdups <- mapToList(SampMap[, c("assay", "colname")])
    logchecks <- any(vapply(lcheckdups, FUN = function(x) {
        as.logical(anyDuplicated(x))
    }, FUN.VALUE = logical(1L)))
    if (logchecks) {
        return("All sample identifiers in the assays must be unique")
    }
    NULL
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
#' @param object A \code{MultiAssayExperiment} object
setMethod("show", "MultiAssayExperiment", function(object) {
    if (.hasOldAPI(object)) {
        stop("MultiAssayExperiment is outdated, please run updateObject()")
    }
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    if (length(o_names) == 0L) {
        o_names <- "none"
    }
    classes <- vapply(experiments(object), class, character(1))
    c_elist <- class(experiments(object))
    c_mp <- class(colData(object))
    c_sm <- class(sampleMap(object))
    cat(sprintf("A %s", o_class),
        "object of", o_len, "listed\n",
        ifelse(o_len == 1L, "experiment", "experiments"),
        "with",
        ifelse(all(o_names == "none"), "no user-defined names",
               ifelse(length(o_names) == 1L, "a user-defined name",
                      "user-defined names")),
        ifelse(length(o_len) == 0L, "or", "and"),
        ifelse(length(o_len) == 0L, "classes.",
               ifelse(length(classes) == 1L,
                      "respective class.", "respective classes.")),
        "\n Containing an ")
    show(experiments(object))
    cat("Features: \n experiments() - obtain the",
        sprintf("%s", c_elist), "instance",
        "\n colData() - the primary/phenotype", sprintf("%s", c_mp),
        "\n sampleMap() - the sample availability", sprintf("%s", c_sm),
        "\n `$`, `[`, `[[` - extract colData columns, subset, or experiment",
        "\n *Format() - convert",
        "into a long or wide", sprintf("%s", c_mp),
        "\n assays() - convert", sprintf("%s", c_elist),
        "to a SimpleList of matrices\n")
})

#' @name MultiAssayExperiment-methods
#' @title Accessing/modifying slot information
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
#'    \item pData: (deprecated) Access the \code{colData} slot
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
#' }
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param object A \code{MultiAssayExperiment} object
#' @param name A column in \code{colData}
#' @param value See details.
#' @param ... Argument not in use
#'
#' @return Accessors: Either a \code{sampleMap}, \code{ExperimentList}, or
#' \code{DataFrame} object
#' @return Setters: A \code{MultiAssayExperiment} object
#'
#' @example inst/scripts/MultiAssayExperiment-methods-Ex.R
#'
#' @aliases experiments sampleMap experiments<- sampleMap<-
NULL

### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

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

#' @exportMethod pData
#' @rdname MultiAssayExperiment-methods
setMethod("pData", "MultiAssayExperiment", function(object) {
    .Defunct("colData")
    colData(object)
})

#' @exportMethod colData
#' @rdname MultiAssayExperiment-methods
setMethod("colData", "MultiAssayExperiment", function(x, ...) {
    getElement(x, "colData")
})

#' @exportMethod metadata
#' @rdname MultiAssayExperiment-methods
setMethod("metadata", "MultiAssayExperiment", function(x)
    getElement(x, "metadata"))

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

setGeneric("sampleMap<-", function(object, value) {
    standardGeneric("sampleMap<-")
})

#' @exportMethod sampleMap<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "DataFrame"),
                function(object, value) {
                    slot(object, "sampleMap") <- value
                    return(object)
                })

setGeneric("experiments<-", function(object, value)
    standardGeneric("experiments<-"))

#' @exportMethod experiments<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("experiments", c("MultiAssayExperiment", "ExperimentList"),
                 function(object, value) {
                     if (!length(value)) {
                         slot(object, "ExperimentList") <- value
                         return(object)
                     }
                     rebliss <- .harmonize(value,
                                           colData(object),
                                           sampleMap(object))
                     BiocGenerics:::replaceSlots(object,
                             ExperimentList = rebliss[["experiments"]],
                             colData = rebliss[["colData"]],
                             sampleMap = rebliss[["sampleMap"]],
                             metadata = metadata(object))
                 })

#' @exportMethod pData<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("pData", c("MultiAssayExperiment", "DataFrame"),
    function(object, value) {
        .Defunct("colData")
        slot(object, "colData") <- value
        return(object)
    })

#' @exportMethod colData<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("colData", c("MultiAssayExperiment", "DataFrame"),
    function(x, value) {
        slot(x, "colData") <- value
        return(x)
    })

.rearrangeMap <- function(sampMap) {
    return(DataFrame(assay = factor(sampMap[["assayname"]]),
                     primary = sampMap[["primary"]],
                     colname = sampMap[["assay"]]))
}

#' @exportMethod metadata<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("metadata", c("MultiAssayExperiment", "ANY"),
                 function(x, ..., value) {
                     slot(x, "metadata") <- value
                     return(x)
                 })

#' @exportMethod $<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("$", "MultiAssayExperiment", function(x, name, value) {
    colData(x)[[name]] <- value
    return(x)
})

#' @exportMethod updateObject
#' @param verbose logical (default FALSE) whether to print extra messages
#' @describeIn MultiAssayExperiment Update old serialized MultiAssayExperiment
#' objects to new API
setMethod("updateObject", "MultiAssayExperiment",
    function(object, ..., verbose = FALSE) {
        if (verbose)
            message("updateObject(object = 'MultiAssayExperiment')")
        oldAPI <- try(object@ExperimentList, silent = TRUE)
        new(class(object),
            ExperimentList = if (is(oldAPI, "try-error"))
                ExperimentList(object@Elist@listData)
            else experiments(object),
            colData = if (is(try(object@colData, silent = TRUE),
                             "try-error"))
                object@pData
            else
                colData(object),
            sampleMap = if (is(oldAPI, "try-error"))
                .rearrangeMap(sampleMap(object))
            else sampleMap(object),
            metadata = metadata(object),
            drops = object@drops)
    })

### ==============================================
### MatchedAssayExperiment class
### ----------------------------------------------

#' An integrative and matched-samples class for experiment data
#'
#' @description
#' This class supports the use of matched samples where an equal number
#' of observations per biological unit are present in all assays.
#'
#' @return A \code{MatchedAssayExperiment} object
#'
#' @exportClass MatchedAssayExperiment
#' @seealso MultiAssayExperiment
setClass("MatchedAssayExperiment", contains="MultiAssayExperiment")

.checkEqualPrimaries <- function(object) {
    listMap <- mapToList(sampleMap(object))
    primaryIDs <- lapply(listMap, function(x) x[["primary"]])
    allIDsEqual <- all(vapply(seq_along(primaryIDs)[-1], function(i, prim) {
        identical(prim[[1L]], prim[[i]])
    }, FUN.VALUE = logical(1L), prim = primaryIDs))
    if (!allIDsEqual)
        "Primary identifiers are not equal across assays"
    else
        NULL
}

.checkPrimaryOrder <- function(object) {
    colPrimary <- rownames(colData(object))
    listMap <- mapToList(sampleMap(object))
    primaryIDs <- lapply(listMap, function(x) x[["primary"]])
    allOrdered <- all(vapply(primaryIDs, function(prim) {
        identical(colPrimary, prim)
    }, logical(1L)))
    if (!allOrdered)
        "colData row identifiers not identical to sampleMap primary column"
    else
        NULL
}

.validMatchedAssayExperiment <- function(object) {
    if (length(object) != 0L) {
    c(.checkEqualPrimaries(object),
      .checkPrimaryOrder(object))
    }
}

S4Vectors::setValidity2("MatchedAssayExperiment", .validMatchedAssayExperiment)

#' @describeIn MatchedAssayExperiment-class Construct a
#' \code{MatchedAssayExperiment} class from \linkS4class{MultiAssayExperiment}
#' inputs.
#'
#' @inheritParams MultiAssayExperiment
#'
#' @examples
#' data("miniACC")
#' acc <- as(miniACC, "MatchedAssayExperiment")
#' acc
#'
#' @aliases MatchedAssayExperiment
#' @export MatchedAssayExperiment
MatchedAssayExperiment <- function(experiments = ExperimentList(),
    colData = S4Vectors::DataFrame(), sampleMap =
        S4Vectors::DataFrame(assay = factor(), primary = character(),
                             colname = character()),
    metadata = NULL, drops = list()) {
    matched <- MultiAssayExperiment(experiments = experiments, colData = colData,
                             sampleMap = sampleMap, metadata = metadata,
                             drops = drops)
    if (!isEmpty(matched))
        matched <- as(matched, "MatchedAssayExperiment")
    return(matched)
}

setAs("MultiAssayExperiment", "MatchedAssayExperiment", function(from) {
    from <- intersectColumns(from)
    if (length(Filter(length, duplicated(from))))
        stop("Resolve replicate columns")
    newMatched <- new("MatchedAssayExperiment", from)
    return(newMatched)
})
