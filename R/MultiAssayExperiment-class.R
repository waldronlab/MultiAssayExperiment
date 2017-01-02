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

#' An integrative MultiAssay class for experiment data
#'
#' @description
#' The \code{MultiAssayExperiment} class can be used to manage results of
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of
#' \code{\linkS4class{SummarizedExperiment}},
#' \code{\linkS4class{ExpressionSet}},
#' \code{matrix}, \code{\link{RangedRaggedAssay}}
#' (inherits from \code{\linkS4class{GRangesList}}), and \code{RangedVcfStack}.
#' Create new \code{MultiAssayExperiment} instances with the eponymous
#' constructor, minimally with the argument \code{\link{ExperimentList}},
#' potentially also with the arguments \code{pData} (see section below) and
#' \code{\link{sampleMap}}.
#'
#' @details
#' The dots (\code{\ldots}) argument allows the user to specify additional
#' arguments in serveral instances. When subsetting a
#' \code{MultiAssayExperiment}, the dots allow for additional
#' arguments to be sent to \link{findOverlaps}. When using the \code{reduce}
#' method, the dots are used to specify arguments for the supplied
#' \code{combine} argument and function.
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
#' MultiAssayExperiment()
#'
#' @exportClass MultiAssayExperiment
#' @include ExperimentList-class.R
#' @seealso SubsetByRow method: \link{getHits}
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
    cat("To access: \n experiments() - to obtain the",
        sprintf("%s", c_elist), "instance",
        "\n pData() - for the primary/phenotype", sprintf("%s", c_mp),
        "\n sampleMap() - for the sample availability", sprintf("%s", c_sm),
        "\n metadata() - for the metadata object of ANY class",
        "\nSee also: subsetByAssay(), subsetByRow(), subsetByColumn()\n")
})


### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

#' Accessor function for the \code{sampleMap} slot of a
#' \code{MultiAssayExperiment} object
#'
#' @param x A \code{MultiAssayExperiment} object
#' @return A \code{DataFrame} object of sample relationships across experiments
#' @example inst/scripts/sampleMap-Ex.R
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))

#' @describeIn MultiAssayExperiment Access sampleMap slot from a
#' MultiAssayExperiment
#' @exportMethod sampleMap
setMethod("sampleMap", "MultiAssayExperiment", function(x)
    getElement(x, "sampleMap"))

#' Accessor function for the \code{ExperimentList} slot of a
#' \code{MultiAssayExperiment} object
#'
#' @param x A \code{MultiAssayExperiment} class object
#' @return A \code{ExperimentList} object of assay data
#' @examples
#' example("MultiAssayExperiment")
#' experiments(myMultiAssayExperiment)
setGeneric("experiments", function(x) standardGeneric("experiments"))

#' @describeIn MultiAssayExperiment Access ExperimentList class from a
#' MultiAssayExperiment
#' @exportMethod experiments
setMethod("experiments", "MultiAssayExperiment", function(x)
    getElement(x, "ExperimentList"))

#' @describeIn MultiAssayExperiment Access pData slot from a
#' MultiAssayExperiment
#' @exportMethod pData
#' @importFrom Biobase pData
setMethod("pData", "MultiAssayExperiment", function(object)
    getElement(object, "pData"))

#' @describeIn MultiAssayExperiment Access metadata slot from a
#' MultiAssayExperiment
#' @exportMethod metadata
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

#' Replace a slot value with a given \code{DataFrame}
#'
#' @param object A \code{MultiAssayExperiment} object
#' @param value A \code{DataFrame} object to replace the existing
#' \code{sampleMap}
#'
#' @examples
#' ## Load example
#' example("MultiAssayExperiment")
#'
#' ## Replacement method for a MultiAssayExperiment sampleMap
#' sampleMap(myMultiAssayExperiment) <- DataFrame()
#'
#' @return A \code{sampleMap} with replacement values
setGeneric("sampleMap<-", function(object, value) {
    standardGeneric("sampleMap<-")
})

#' @exportMethod sampleMap<-
#' @describeIn MultiAssayExperiment value: A \code{DataFrame} sampleMap
#' representation
#' @param value A \code{DataFrame} or \code{ExperimentList} object to replace
#' the existing \code{sampleMap}, \code{ExperimentList}, or \code{pData} slot
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "DataFrame"),
                function(object, value) {
                    slot(object, "sampleMap") <- value
                    return(object)
                })

#' Replace an \code{ExperimentList} slot value with a given
#' \code{ExperimentList} class object
#'
#' @param object A \code{MultiAssayExperiment} class object
#' @param value An \code{ExperimentList} object to replace the existing
#' \code{ExperimentList} slot
#'
#' @examples
#' ## Load a MultiAssayExperiment
#' example("MultiAssayExperiment")
#'
#' ## Replace with an empty ExperimentList
#' experiments(myMultiAssayExperiment) <- ExperimentList()
#'
#' @return A \code{ExperimentList} class object
setGeneric("experiments<-", function(object, value)
    standardGeneric("experiments<-"))

#' @exportMethod experiments<-
#' @describeIn MultiAssayExperiment value: An \code{ExperimentList}
#' representation
setReplaceMethod("experiments", c("MultiAssayExperiment", "ExperimentList"),
                function(object, value) {
                    slot(object, "ExperimentList") <- value
                    return(object)
                })

#' @exportMethod pData<-
#' @describeIn MultiAssayExperiment value: A \code{DataFrame} of specimen data
#' @importFrom Biobase pData<-
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
#' @describeIn MultiAssayExperiment value: Data of type "ANY"
#' @importFrom S4Vectors metadata<-
setReplaceMethod("metadata", c("MultiAssayExperiment", "ANY"),
                 function(x, ..., value) {
                     slot(x, "metadata") <- value
                     return(x)
                 })

#' @exportMethod $<-
#' @describeIn MultiAssayExperiment value: DataFrame column
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
              return(object)
          })
