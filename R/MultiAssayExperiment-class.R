### ==============================================
### MultiAssayExperiment class
### ----------------------------------------------

#' An integrative MultiAssay class for experiment data
#' 
#' @description 
#' The \code{MultiAssayExperiment} class can be used to manage results of 
#' diverse assays on a collection of specimen. Currently,  the class can handle
#' assays that are organized instances of \code{RangedSummarizedExperiment}, 
#' \code{ExpressionSet}, \code{matrix}, \code{RangedRaggedAssay} 
#' (inherits from \code{GRangesList}), and \code{RangedVcfStack}. Create new
#' \code{MultiAssayExperiment} instances with the eponymous constructor, 
#' minimally with the argument \code{\linkS4class{Elist}}, potentially also
#' with the arguments \code{pData} and \code{sampleMap}
#' 
#' @slot Elist A \code{\linkS4class{Elist}} class object for each assay dataset
#' @slot pData A \code{DataFrame} of all clinical data available across
#' experiments
#' @slot sampleMap A \code{DataFrame} of translatable identifiers of samples
#' and participants
#' @slot metadata Additional data describing the
#' \code{\link{MultiAssayExperiment}} object
#' @slot drops A metadata \code{list} of dropped information.
#' @exportClass MultiAssayExperiment
#' @include Elist-class.R
setClass("MultiAssayExperiment",
         slots = list(
           Elist = "Elist",
           pData = "DataFrame", 
           sampleMap = "DataFrame",
           metadata = "ANY",
           drops = "list"
         )
)

##
## Validity ---------------------------------
##

.uniqueSortIdentical <- function(charvec1, charvec2) {
  listInput <- list(charvec1, charvec2)
  listInput <- lapply(listInput, function(x) sort(unique(x)))
  return(identical(listInput[[1]], listInput[[2]]))
}

.allIn <- function(charvec1, charvec2) {
  return(all(charvec2 %in% charvec1))
}

## sampleMap is a DataFrame with unique sampleNames across assay
.checkSampleMapNames <- function(object) {
  errors <- character()
  if (!(.allIn(rownames(pData(object)),
               slot(sampleMap(object)[, "master"], "values")))) {
    msg <- "All samples in the sampleMap must be in the pData"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0L) 
    NULL else errors
}

.uniqueNamesInAssays <- function(object) {
  errors <- character()
  SampMap <- sampleMap(object)
  lcheckdups <- S4Vectors::split(SampMap[["assay"]], SampMap[, "assayname"])
  logchecks <- any(vapply(lcheckdups,
                          function(x) any(duplicated(x)), logical(1)))
  if (logchecks) {
    msg <- "All sample identifiers in the assays must be unique"
    errors <- c(errors, msg)
  }
  if (length(errors) == 0L) 
    NULL else errors
}

## Experiment list must be the same length as the unique sampleMap assaynames
.checkElist2 <- function(object) {
  errors <- character()
  assaynames <- unique(sampleMap(object)[, "assayname"])
  if (length(Elist(object)) != length(assaynames)) {
    msg <- "Elist must be the same length as the sampleMap assaynames"
    errors <- c(errors, msg)
  }
  if (!all(assaynames %in% names(object))) {
    msg <- paste0("Experiment/Assay names in both the ",
                  "Elist and the sampleMap must match")
    errors <- c(errors, msg)
  }
  if (length(errors) == 0L) 
    NULL else errors
}

## All sample names in the Elist must be in the sampleMap
.checkSampleNames <- function(object) {
  if (!identical(sort(unname(unlist(colnames(object)))),
                 sort(sampleMap(object)[, "assay"]))) {
    return("samples in the Elist are not the same as samples in the sampleMap")
  }
  NULL
}

## All names must match between Elist and sampleMap
.checkNames <- function(object) {
  if (!all(
    names(Elist(object)) %in% unique(sampleMap(object)[, "assayname"]))) {
    return("Experiment names must match in both Elist and sampleMap")
  }
  NULL
}

.validMultiAssayExperiment <- function(object) {
  if (length(Elist(object)) != 0L) {
    c(.checkNames(object),
      .checkElist2(object),
      .checkSampleMapNames(object), 
      .uniqueNamesInAssays(object),
      .checkSampleNames(object)
     )
  }
}

S4Vectors::setValidity2("MultiAssayExperiment", .validMultiAssayExperiment)


#' @exportMethod show
#' @describeIn MultiAssayExperiment Show method for a
#' \code{MultiAssayExperiment}
#' @param object A \code{MultiAssayExperiment} class object
setMethod("show", "MultiAssayExperiment", function(object) {
  o_class <- class(object)
  o_len <- length(object)
  o_names <- names(object)
  if (length(o_names) == 0L) {
    o_names <- "none"
  }
  classes <- vapply(Elist(object), class, character(1))
  c_elist <- class(Elist(object))
  c_mp <- class(pData(object))
  c_sm <- class(sampleMap(object))
  cat(sprintf('A "%s"', o_class),
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
  show(Elist(object))
  cat("To access slots use: \n Elist() - to obtain the",
      sprintf('"%s"', c_elist), 
      "of experiment instances",
      "\n pData() - for the phenotype", sprintf('"%s"', c_mp),
      "\n sampleMap() - for the sample availability", sprintf('"%s"', c_sm),
      "\n metadata() - for the metadata object of 'ANY' class",
      "\nSee also: subsetByAssay(), subsetByFeature(), subsetBySample()\n")
})


### - - - - - - - - - - - - - - - - - - - - - - -
### Accessor methods
###

#' Accessor function for the \code{sampleMap} slot of a
#' \code{MultiAssayExperiment} object
#'
#' @param x A \code{MultiAssayExperiment} object
#' @return A \code{DataFrame} object of sample relationships across experiments
setGeneric("sampleMap", function(x) standardGeneric("sampleMap"))
#' @describeIn MultiAssayExperiment Access sampleMap slot from
#' MultiAssayExperiment
#' @exportMethod sampleMap
setMethod("sampleMap", "MultiAssayExperiment", function(x)
  getElement(x, "sampleMap"))

#' @describeIn MultiAssayExperiment Access Elist class from
#' MultiAssayExperiment
#' @exportMethod Elist
setMethod("Elist", "MultiAssayExperiment", function(x)
  getElement(x, "Elist"))

#' @describeIn MultiAssayExperiment Access pData slot from
#' MultiAssayExperiment
#' @exportMethod pData
#' @importFrom Biobase pData
setMethod("pData", "MultiAssayExperiment", function(object)
  getElement(object, "pData"))

#' @describeIn MultiAssayExperiment Access metadata slot from
#' MultiAssayExperiment
#' @exportMethod metadata
setMethod("metadata", "MultiAssayExperiment", function(x)
  getElement(x, "metadata"))

### - - - - - - - - - - - - - - - - - - - - - - - -
### Getters
###

#' @exportMethod length
#' @describeIn MultiAssayExperiment Get the length of Elist 
setMethod("length", "MultiAssayExperiment", function(x)
  length(getElement(x, "Elist"))
)

#' @exportMethod names
#' @describeIn MultiAssayExperiment Get the names of the Elist
setMethod("names", "MultiAssayExperiment", function(x)
  names(getElement(x, "Elist"))
)

#' Replace a slot value with a given \code{DataFrame}
#' 
#' @param x A \code{MultiAssayExperiment} object
#' @param value A \code{DataFrame} object to replace the existing
#' \code{sampleMap}
setGeneric("sampleMap<-", function(x, value) standardGeneric("sampleMap<-"))
#' @exportMethod sampleMap<-
#' @describeIn MultiAssayExperiment value: A \code{DataFrame} sampleMap 
#' representation
#' @param value A \code{DataFrame} or \code{Elist} object to replace the existing
#' \code{sampleMap} or an \code{Elist} slot, respectively
setReplaceMethod("sampleMap", c("MultiAssayExperiment", "DataFrame"),
                 function(x, value) {
                   slot(x, "sampleMap") <- value
                   return(x)
                 })
#' Replace an \code{Elist} slot value with a given \code{\linkS4class{Elist}}
#' class object
#'
#' @param x A \code{MultiAssayExperiment} class object
#' @param value An \code{\linkS4class{Elist}} object to replace the existing
#' \code{\linkS4class{Elist}} slot
setGeneric("Elist<-", function(x, value) standardGeneric("Elist<-"))
#' @exportMethod Elist<-
#' @describeIn MultiAssayExperiment value: An \linkS4class{Elist} 
#' representation
setReplaceMethod("Elist", c("MultiAssayExperiment", "Elist"),
                 function(x, value) {
                   slot(x, "Elist") <- value
                   return(x)
                 })
