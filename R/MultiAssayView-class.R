### ==============================================
### MultiAssayView class
### ==============================================

#' An identifier class used for staging a subset operation
#'
#' @slot query Any class indicator needed to subset
#' @slot keeps A \code{list} indicating valid matches in each assay
#' @slot drops A \code{list} of excluded information due to subsetting
#' @slot type A \code{character} vector indicating method used to search
#' @exportClass MultiAssayView

.MultiAssayView <- setClass("MultiAssayView",
         representation(
             subject = "environment",
             sampleMap="DataFrame",
             features="list")
         )

MultiAssayView <- function(subject) {
    .MultiAssayView(subject=as.environment(list(subject=subject)),
                    sampleMap=sampleMap(subject))
}

setMethod("sampleMap", "MultiAssayView", function(x) {
    slot(x, "sampleMap")
})

## FIXME
setGeneric("sampleMap<-", function(x, value) standardGeneric("sampleMap<-"))

setReplaceMethod("sampleMap", c("MultiAssayView", "DataFrame"),
    function(x, value)
{
    slot(x, "sampleMap") <- value
    x
})

setMethod("colnames", "MultiAssayView", function(x) {
    map <- sampleMap(x)
    split(map$assay, map$assayname)
})

setMethod("rownames", "MultiAssayView", function(x) {
    ## FIXME
})

setMethod("[", c("MultiAssayView", "missing", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    ## FIXME: what are we subsetting on?
    stopifnot(all(j %in% sampleMap(x)$master))
    map <- sampleMap(x)
    sampleMap(x) <- map[match(j, map$master),, drop=FALSE]
    x
})

setMethod("show", "MultiAssayView", function(object) {
    map <- sampleMap(object)
    assays <- unique(map$assayname)
    cat("A ", class(object), " object",
        "\nExperiments:\n  ",
        paste(assays, collapse="\n  "),
        "\nSamples: ", length(unique(map$master)),
        "\n", sep="")
})
