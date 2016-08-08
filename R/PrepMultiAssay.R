#' Prepare a \code{MultiAssayExperiment} instance
#'
#' The purpose of this helper function is to faciltate the creation of a
#' \code{\link{MultiAssayExperiment}} object by detecting any inconsistencies
#' with all types of names in either the \code{\link{ExperimentList}},
#' the \code{pData}, or \code{\link{sampleMap}}.
#'
#' @section Checks:
#' The \code{PrepMultiAssay} function checks that all columns in the sampleMap
#' are \code{character}.
#'
#' It checks that all names and lengths match in both the
#' \code{\link{ExperimentList}} and in the unique assaynames of the
#' \code{\link{sampleMap}}.
#'
#' If \code{\link{ExperimentList}} names and assaynames only differ by case
#' and are not #' duplicated, the function will standardize all names to
#' lowercase.
#'
#' If names cannot be matched between the assay column of the
#' \code{\link{sampleMap}} and the colnames of the \code{ExperimentList}, those
#' unmatched will be dropped and found in the "drops" element of the
#' resulting \code{list}.
#'
#' Names in the "primary" column of the \code{\link{sampleMap}}, will be
#' matched to those in the \code{pData}. Unmatched "primary" column rows will
#' be dropped from the \code{\link{sampleMap}}. Suggestions for name fixes in
#' either the \code{\link{ExperimentList}} or colnames will be made when
#' necessary.
#'
#' @param ExperimentList A \code{list} of all combined experiments
#' @param pData A \linkS4class{DataFrame} of the phenotype
#' data for all participants
#' @param sampleMap A \linkS4class{DataFrame} of sample identifiers, assay
#' samples, and assay names
#' @return A \code{list} containing all the essential components of a
#' \code{\link{MultiAssayExperiment}} as well as a "drops" element that
#' indicates non-matched names.
#'
#' @examples
#' ## Run example
#' example("MultiAssayExperiment")
#'
#' ## Check if there are any inconsistencies within the different names
#' preparedMAE <- PrepMultiAssay(ExpList, pDat, mySampleMap)
#'
#' ## Results in a list of components for the MultiAssayExperiment constructor
#' ## function
#' MultiAssayExperiment(preparedMAE$ExperimentList, preparedMAE$pData,
#' preparedMAE$sampleMap)
#'
#' @export PrepMultiAssay
PrepMultiAssay <- function(ExperimentList, pData, sampleMap) {
    drops <- list()
    ExperimentList <- ExperimentList(ExperimentList)
    if (any(vapply(sampleMap, FUN = function(col) {
        !is.character(col)
    }, FUN.VALUE = logical(1L)))) {
        sampleMap[] <- lapply(sampleMap, as.character)
    }
    if (!is(sampleMap, "DataFrame"))
        sampleMap <- S4Vectors::DataFrame(sampleMap)
    if (is.null(names(ExperimentList)))
        stop("ExperimentList does not have names, assign names")
    assaynames <- unique(sampleMap[["assay"]])
    if (length(names(ExperimentList)) != length(assaynames)) {
        warning("Lengths of names in the ExperimentList and sampleMap",
                " are not equal")
    } else if (any(!(assaynames %in% names(ExperimentList)))) {
        message("Names in the ExperimentList do not match sampleMap assaynames",
                "\nstandardizing will be attempted...")
        nameErr <- TRUE
        if (identical(tolower(assaynames), tolower(names(ExperimentList))) &&
            !as.logical(anyDuplicated(tolower(assaynames),
                                      tolower(names(ExperimentList))))) {
            message(" - names set to lowercase")
            sampleMap[["assay"]] <- tolower(sampleMap[["assay"]])
            names(ExperimentList) <- tolower(names(ExperimentList))
            nameErr <- FALSE
        } else {
            warning("ExperimentList and sampleMap assaynames are not equal")
        }
    }
    primaries <- sampleMap[["primary"]]
    notFounds <- primaries %in% rownames(pData)
    if (!all(notFounds)) {
        message("Not all names in the primary column of the sampleMap",
                "\n  could be matched to the pData rownames; see $drops")
        notF <- sampleMap[!notFounds, ]
        drops <- list(sampleMap_rows = notF)
        sampleMap <- sampleMap[notFounds,]
        print(notF)
        if (length(unique(sampleMap[["assay"]])) != length(ExperimentList)) {
            stop("Some assays could not be matched,",
                 " check primary and pData names")
        }
    }
    if (exists("nameErr") && nameErr) {
        stop("Fix ExperimentList and sampleMap assaynames before checking for ",
             "matchable columns")
    }
    listMap <- mapToList(sampleMap)
    cols <- colnames(ExperimentList)
    listMap <- listMap[names(ExperimentList)]
    allThere <- mapply(function(colnams, sampmap) {
        colnams %in% sampmap[, "colname"]
    }, colnams = cols, sampmap = listMap, SIMPLIFY = FALSE)
    whichNotAll <- vapply(allThere, FUN = function(logicalVector) {
        !(all(logicalVector))
    }, FUN.VALUE = logical(1L))
    if (any(whichNotAll)) {
        message("Not all colnames in the ExperimentList are found in the \n",
                "sampleMap, dropping samples from ExperimentList...")
        ExperimentList <- ExperimentList(mapply(function(x, y) {
            x[, y, drop = FALSE]
        }, x = ExperimentList, y = allThere, SIMPLIFY = FALSE))
        coldrops <- mapply(function(listColnames, logicalList) {
            listColnames[!logicalList]
        }, listColnames = cols, logicalList = allThere, SIMPLIFY = FALSE)
        print(Biobase::selectSome(coldrops))
        drops <- c(drops, columns = coldrops)
    }
    return(list(ExperimentList = ExperimentList, pData = pData,
                sampleMap = sampleMap, drops = drops))
}
