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
#' \code{matrix}, \code{\link{RangedRaggedAssay}}
#' (inherits from \code{\linkS4class{GRangesList}}), and \code{RangedVcfStack}.
#' Create new \code{MultiAssayExperiment} instances with the homonymous
#' constructor, minimally with the argument \code{\link{ExperimentList}},
#' potentially also with the arguments \code{pData} (see section below) and
#' \code{\link{sampleMap}}.
#'
#' @details
#' The dots (\code{\ldots}) argument allows the user to specify additional
#' arguments in several instances. When subsetting (\strong{[}) a
#' \code{MultiAssayExperiment}, the dots allow for additional
#' arguments to be sent to \link{findOverlaps}. When using the \strong{reduce}
#' method, the dots are used to specify arguments for the supplied
#' \code{combine} argument and function. When using the \strong{assay} method,
#' additional arguments may be passed to the \code{RangedRaggedAssay} method.
#' See the link for more information:
#' \link{assay,RangedRaggedAssay,missing-method}. When using \strong{c} method
#' to add experiments to a \code{MultiAssayExperiment}, the dots allow extra
#' data classes compatible with the MultiAssayExperiment API. See: \link{API}
#'
#' @section pData:
#' The \code{pData} slot is a collection of primary specimen data valid across
#' all experiments. This slot is strictly of class
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
#' @slot pData A \code{DataFrame} of all clinical/specimen data available
#' across experiments
#' @slot sampleMap A \code{DataFrame} of translatable identifiers
#' of samples and participants
#' @slot metadata Additional data describing the
#' \code{MultiAssayExperiment} object
#' @slot drops A metadata \code{list} of dropped information
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
#' # Columns (j) Rows in pData
#' myMultiAssayExperiment[, rownames(pData(myMultiAssayExperiment))[3:2],  ]
#'
#' # Assays (k)
#' myMultiAssayExperiment[, , "Affy"]
#'
#' ## Complete cases (returns logical vector)
#' completes <- complete.cases(myMultiAssayExperiment)
#' compMAE <- myMultiAssayExperiment[, completes, ]
#' compMAE
#' pData(compMAE)
#'
#' @exportClass MultiAssayExperiment
#' @seealso \link{MultiAssayExperiment-methods} for slot modifying methods
#' @include ExperimentList-class.R
setClass("MultiAssayExperiment",
         slots = list(
           ExperimentList = "ExperimentList",
           pData = "DataFrame",
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

## PDATA
## 2.i. See setClass above where pData = "DataFrame"

## SAMPLEMAP
## 3.i. all values in the sampleMap "primary" column must be found in the
## rownames of pData
.checkSampleMapNames <- function(object) {
    errors <- character()
    if (!(.allIn(
        rownames(pData(object)),
        sampleMap(object)[["primary"]]
    ))) {
        msg <- "All samples in the sampleMap must be in the pData"
        errors <- c(errors, msg)
    }
    if (length(errors) == 0L)
        NULL else errors
}

## 3.ii. Within rows of "sampleMap" corresponding to a single value in the
## "assayname" column, there can be no duplicated values in the "assay" column
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
          .checkSampleMapNames(object),
          .uniqueNamesInAssays(object),
          .checkSampleNames(object)
        )
    }
}

S4Vectors::setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)

.hasOldAPI <- function(object) {
    isTRUE(.hasSlot(object, "Elist"))
}

#' @exportMethod show
#' @describeIn MultiAssayExperiment Show method for a
#' \code{MultiAssayExperiment}
#' @param object A \code{MultiAssayExperiment} class object
setMethod("show", "MultiAssayExperiment", function(object) {
    if (.hasOldAPI(object)) {
        object <- updateObject(object)
        warning("MultiAssayExperiment is outdated, please run updateObject()")
    }
    o_class <- class(object)
    o_len <- length(object)
    o_names <- names(object)
    if (length(o_names) == 0L) {
        o_names <- "none"
    }
    classes <- vapply(experiments(object), class, character(1))
    c_elist <- class(experiments(object))
    c_mp <- class(pData(object))
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
        "\n pData() - the primary/phenotype", sprintf("%s", c_mp),
        "\n sampleMap() - the sample availability", sprintf("%s", c_sm),
        "\n `$`, `[`, `[[` - extract pData columns, subset, or experiment",
        "\n reduce() - select complete cases, order columns, disjoin ranges",
        "\n rearrange() - convert", sprintf("%s", c_elist),
        "into a long or wide", sprintf("%s", c_mp),
        "\n assay() - convert", sprintf("%s", c_elist),
        "to a list of rectangular matrices\n")
})

#' @name MultiAssayExperiment-methods
#' @title Accessing/modifying slot information
#'
#' @description A set of accessor and setter generic functions to extract
#' either the \code{sampleMap}, the \code{\link{ExperimentList}}, \code{pData},
#' or \code{metadata} slots of a \code{\link{MultiAssayExperiment}} object
#'
#' @section Accessors:
#' Eponymous names for accessing \code{MultiAssayExperiment} slots with the
#' exception of the \link{ExperimentList} accessor named \code{experiments}.
#' \itemize{
#'    \item experiments: Access the \link{ExperimentList} slot
#'    \item [[: Access the \link{ExperimentList} slot
#'    \item $: Access a column in \code{pData}
#' }
#'
#' @section Setters:
#' Setter method values (i.e., '\code{function(x) <- value}'):
#' \itemize{
#'     \item experiments<-: An \code{\link{ExperimentList}} object
#'     containing experiment data of supported classes
#'     \item sampleMap<-: A \code{\link{DataFrame}} object relating
#'     samples to biological units and assays
#'     \item pData<-: A \code{\link{DataFrame}} object describing the
#'     biological units
#'     \item metadata<-: A \code{list} object of metadata
#'     \item [[<-: Equivalent to the '\code{experiments<-}' setter method for
#'     convenience
#'     \item $<-: A vector to replace the indicated column in \code{pData}
#' }
#'
#' @param x A \code{MultiAssayExperiment} object
#' @param object A \code{MultiAssayExperiment} object
#' @param name A column in \code{pData}
#' @param value See details.
#' @param i A \code{numeric} or \code{character} vector of length 1
#' @param j Argument not in use
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
#'
#' @importFrom Biobase pData
setMethod("pData", "MultiAssayExperiment", function(object)
    getElement(object, "pData"))

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
    length(getElement(x, "ExperimentList"))
)

#' @exportMethod names
#' @describeIn MultiAssayExperiment Get the names of the ExperimentList
setMethod("names", "MultiAssayExperiment", function(x)
    names(getElement(x, "ExperimentList"))
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
                    slot(object, "ExperimentList") <- value
                    return(object)
                })

#' @exportMethod pData<-
#' @importFrom Biobase pData<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("pData", c("MultiAssayExperiment", "DataFrame"),
                function(object, value) {
                 slot(object, "pData") <- value
                 return(object)
                })

.rearrangeMap <- function(sampMap) {
    return(DataFrame(assay = sampMap[["assayname"]],
                     primary = sampMap[["primary"]],
                     colname = sampMap[["assay"]]))
}

#' @exportMethod metadata<-
#' @importFrom S4Vectors metadata<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("metadata", c("MultiAssayExperiment", "ANY"),
                 function(x, ..., value) {
                     slot(x, "metadata") <- value
                     return(x)
                 })

#' @exportMethod $<-
#' @rdname MultiAssayExperiment-methods
setReplaceMethod("$", "MultiAssayExperiment", function(x, name, value) {
    pData(x)[[name]] <- value
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
              if (is(try(object@ExperimentList, silent = TRUE), "try-error")) {
                  object <- new(class(object),
                                ExperimentList = ExperimentList(
                                    object@Elist@listData),
                                pData = pData(object),
                                sampleMap = .rearrangeMap(sampleMap(object)),
                                metadata = metadata(object),
                                drops = object@drops)
              }
              classes <- vapply(experiments(object), class, character(1L))
              if (any(classes %in% "RangedRaggedAssay")) {
                  rraIdx <- which(classes == "RangedRaggedAssay")
                  for (i in rraIdx) {
                      object[[i]] <- as(object[[i]], "RaggedExperiment")
                  }
              }
              return(object)
          })
