#' Prepare a \code{MultiAssayExperiment} instance
#' 
#' The purpose of this helper function is to faciltate the creation of a
#' \code{\link{MultiAssayExperiment}} object by detecting any inconsistencies
#' with all types of names in either the \code{\link{Elist}}, the \code{pData},
#' or \code{\link{sampleMap}}.
#' 
#' @section Checks:
#' The \code{PrepMultiAssay} function checks that all columns in the sampleMap
#' are \code{character}.
#' 
#' It checks that all names and lengths match in both the \code{\link{Elist}}
#' and in the unique assaynames of the \code{\link{sampleMap}}.
#' 
#' If \code{\link{Elist}} names and assaynames only differ by case and are not
#' duplicated, the function will standardize all names to lowercase.
#' 
#' If names cannot be matched between the assay column of the
#' \code{\link{sampleMap}} and the colnames of the \code{Elist}, those
#' unmatched will be dropped and found in the "drops" element of the
#' resulting \code{list}.
#' 
#' Names in the "primary" column of the \code{\link{sampleMap}}, will be
#' matched to those in the \code{pData}. Unmatched "primary" column rows will
#' be dropped from the \code{\link{sampleMap}}. Suggestions for name fixes in
#' either the \code{\link{Elist}} or colnames will be made when necessary.
#' 
#' @param Elist A \code{list} of all combined experiments
#' @param pData A \linkS4class{DataFrame} of the phenotype
#' data for all participants
#' @param sampleMap A \linkS4class{DataFrame} of sample identifiers, assay
#' samples, and assay names
#' @return A \code{list} containing all the essential components of a
#' \code{\link{MultiAssayExperiment}} as well as a "drops" element that
#' indicates non-matched names.
#' 
#' @export PrepMultiAssay
PrepMultiAssay <- function(Elist, pData, sampleMap) {
  drops <- list()
  Elist <- Elist(Elist)
  if (any(vapply(sampleMap, FUN = function(col) {
    !is.character(col)
  }, FUN.VALUE = logical(1L)))) {
    sampleMap[] <- lapply(sampleMap, as.character)
  }
  if (!is(sampleMap, "DataFrame")) {
    sampleMap <- S4Vectors::DataFrame(sampleMap)
  }
  if (is.null(names(Elist))) {
    warning("Elist does not have names, assign names")
  }
  assaynames <- unique(sampleMap[, "assayname"])
  if (length(names(Elist)) != length(assaynames)) {
    warning("Lengths of names in the Elist and sampleMap are not equal")
  } else if (any(!(assaynames %in% names(Elist)))) {
    message("Names in the Elist do not match sampleMap assaynames",
            "\nstandardizing will be attempted...")
    nameErr <- TRUE
    if (identical(tolower(assaynames), tolower(names(Elist))) &&
        !as.logical(anyDuplicated(tolower(assaynames),
                                  tolower(names(Elist))))) {
      message(" - names set to lowercase")
      sampleMap[, "assayname"] <- tolower(sampleMap[, "assayname"])
      names(Elist) <- tolower(names(Elist))
      nameErr <- FALSE
    } else {
      warning("Elist and sampleMap assaynames are not equal")
    }
  }
  primaries <- sampleMap[, "primary"]
  notFounds <- primaries %in% rownames(pData)
  if (!all(notFounds)) {
    message("Not all names in the primary column of the sampleMap",
            "\ncould be matched to the pData rownames; see $drops")
    notF <- setdiff(primaries, rownames(pData))
    drops <- list(rows = notF)
    sampleMap <- sampleMap[notFounds,]
    print(Biobase::selectSome(notF))
    if (length(unique(sampleMap[, "assayname"])) != length(Elist)) {
      stop("Some assays could not be matched, check primary and pData names")
    }
  }
  if (exists("nameErr") && nameErr) {
    stop("Fix Elist and sampleMap assaynames before checking for ",
         "matchable columns")
  }
  listMap <- mapToList(sampleMap)
  cols <- colnames(Elist)
  listMap <- listMap[names(Elist)]
  allThere <- mapply(function(colnams, sampmap) {
    colnams %in% sampmap[, "assay"]
  }, colnams = cols, sampmap = listMap)
  whichNotAll <- vapply(allThere, FUN = function(g) {
    !(all(g))
  }, FUN.VALUE = logical(1L))
  if (any(whichNotAll)) {
    Elist <- Elist(mapply(function(x, y) {
      x[, y, drop = FALSE]
    }, x = Elist, y = allThere))
    coldrops <- mapply(function(a, b) {a[!b]}, a = cols, b = allThere)
    drops <- c(drops, columns = coldrops)
  }
  return(list(Elist = Elist, pData = pData, sampleMap = sampleMap,
              drops = drops))
}