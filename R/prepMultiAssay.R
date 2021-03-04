#' Prepare a \code{MultiAssayExperiment} instance
#'
#' The purpose of this helper function is to faciltate the creation of a
#' \code{\link{MultiAssayExperiment}} object by detecting any inconsistencies
#' with all types of names in either the \code{\link{ExperimentList}},
#' the \code{colData}, or \code{\link{sampleMap}}.
#'
#' @section Checks:
#' The \code{prepMultiAssay} function checks that all columns in the sampleMap
#' are \code{character}.
#'
#' It checks that all names and lengths match in both the
#' \code{\link{ExperimentList}} and in the unique assay names of the
#' \code{\link{sampleMap}}.
#'
#' If \code{\link{ExperimentList}} names and assay names only differ by case
#' and are not duplicated, the function will standardize all names to
#' lowercase.
#'
#' If names cannot be matched between the colname column of the
#' \code{\link{sampleMap}} and the colnames of the \code{ExperimentList}, those
#' unmatched will be dropped and found in the "drops" element of the
#' resulting \code{list}.
#'
#' Names in the "primary" column of the \code{\link{sampleMap}}, will be
#' matched to those in the \code{colData}. Unmatched "primary" column rows will
#' be dropped from the \code{\link{sampleMap}}. Suggestions for name fixes in
#' either the \code{\link{ExperimentList}} or colnames will be made when
#' necessary.
#'
#' @param ExperimentList A \code{list} of all combined experiments
#' @param colData A \linkS4class{DataFrame} of the phenotype
#' data for all participants
#' @param sampleMap A \linkS4class{DataFrame} of sample identifiers, assay
#' samples, and assay names
#' @param ... Optional arguments for the \code{MultiAssayExperiment} constructor
#' function such as \code{metadata} and \code{drops}.
#' @return A \code{list} containing all the essential components of a
#' \code{\link{MultiAssayExperiment}} as well as a "drops" metadata element that
#' indicates non-matched names. The names of the resulting list correspond to
#' the arguments of the \code{MultiAssayExperiment} constructor function.
#'
#' @examples
#' ## Run example
#' example("MultiAssayExperiment")
#'
#' ## Check if there are any inconsistencies within the different names
#' preparedMAE <- prepMultiAssay(ExpList, colDat, sampMap)
#'
#' ## Results in a list of components for the MultiAssayExperiment constructor
#' ## function
#' MultiAssayExperiment(preparedMAE$experiments, preparedMAE$colData,
#' preparedMAE$sampleMap)
#'
#' ## Alternatively, use the do.call function
#' do.call(MultiAssayExperiment, preparedMAE)
#'
#' @export prepMultiAssay
prepMultiAssay <- function(ExperimentList, colData, sampleMap, ...) {
    if (!is(sampleMap, "DataFrame"))
        sampleMap <- S4Vectors::DataFrame(sampleMap)
    if (is.null(names(ExperimentList)))
        stop("ExperimentList does not have names, assign names")

    dotargs <- list(...)
    drops <- if (!is.null(dotargs[["drops"]])) dotargs[["drops"]] else list()
    metadata <- if (!is.null(dotargs[["metadata"]]))
        dotargs[["metadata"]] else list()
    ExperimentList <- ExperimentList(ExperimentList)

    charTypes <- list(
        assay = is.factor,
        primary = is.character,
        colname = is.character)

    correctType <- vapply(seq_along(charTypes),
        function(i, sampMap, char) {
            char[[i]](sampMap[, i])
            }, logical(1L), sampMap = sampleMap, char = charTypes)

    if (any(!correctType)) {
        coerceList <- list(assay = as.factor, primary = as.character,
            colname = as.character)
        fxIdx <- which(!correctType)
        for (i in fxIdx) {
            sampleMap[, i] <- coerceList[[i]](sampleMap[[i]])
        }
    }

    assays <- sort(levels(sampleMap[["assay"]]))
    expNames <- sort(names(ExperimentList))
    if (length(expNames) != length(assays)) {
        warning("\nLengths of names in the ExperimentList and sampleMap\n",
                " are not equal")
    } else if (any(!(assays %in% expNames))) {
        message("\nNames in the ExperimentList do not match sampleMap assay",
                "\nstandardizing will be attempted...")
        nameErr <- TRUE
        if (identical(tolower(assays), tolower(expNames)) &&
            !as.logical(anyDuplicated(tolower(assays), tolower(expNames)))) {
            message(" - names set to lowercase")
            sampleMap[["assay"]] <- tolower(sampleMap[["assay"]])
            names(ExperimentList) <- tolower(names(ExperimentList))
            nameErr <- FALSE
        } else {
            warning("\nExperimentList and sampleMap assay names are not equal")
        }
    }
    primaries <- sampleMap[["primary"]]
    notFounds <- primaries %in% rownames(colData)
    if (!all(notFounds)) {
        message("\nNot all names in the primary column of the sampleMap",
                "\n could be matched to the colData rownames; see $drops")
        notF <- sampleMap[!notFounds, ]
        drops <- list(sampleMap = notF)
        sampleMap <- sampleMap[notFounds, ]
        sampleMap[, 1] <- factor(sampleMap[, 1])
        print(notF)
        if (length(unique(sampleMap[["assay"]])) != length(ExperimentList)) {
            stop("Some assay names could not be matched,",
                 " check primary and colData names")
        }
    }
    if (exists("nameErr") && nameErr) {
        stop("Fix ExperimentList and sampleMap assay names before checking ",
             "for matchable columns")
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
        message("\nNot all colnames in the ExperimentList are found in the \n",
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
    return(list(experiments = ExperimentList, colData = colData,
        sampleMap = sampleMap, metadata = list(metadata, drops = drops)))
}
