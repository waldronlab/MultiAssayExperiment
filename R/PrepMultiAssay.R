#' Prepare a MultiAssayExperiment instance
#' 
#' The purpose of this helper function is to faciltate the creation of a
#' MultiAssayExperiment object by detecting any inconsistencies with all types
#' of names in either the Elist or the pData.
#' 
#' @param Elist A \code{list} of all combined experiments
#' @param pData A \code{\link[S4Vectors]{DataFrame-class}} of the phenotype
#' data for all participants
#' @param sampleMap A \code{DataFrame} of sample identifiers, assay samples,
#' and assay names
#' @return Messages and suggestions for creating a MultiAssayExperiment
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