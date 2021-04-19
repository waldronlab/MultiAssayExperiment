#' @param listmap A named \code{list} object containing \code{DataFrame}s
#' with "primary" and "colname" columns
#'
#' @return A \linkS4class{DataFrame} class object of names
#' @describeIn mapToList The inverse of the listToMap operation
#' @export listToMap
listToMap <- function(listmap) {
    if (is.null(names(listmap)))
        stop("'listmap' must be a named list")

    elementClass <- unique(vapply(listmap, class, character(1L)))

    if (length(elementClass) != 1L)
        stop("All 'listmap' elements must be of the same class")

    if (!(is(listmap[[1]], "DataFrame") || is.data.frame(listmap[[1]])))
        stop("'listmap' elements are not all 'DataFrame' or 'data.frame'")

    if (elementClass == "data.frame")
        listmap <- lapply(listmap, DataFrame)

    listmap <- lapply(listmap, function(lmap) {
        if (isEmpty(lmap))
            DataFrame(primary = NA_character_, colname = NA_character_)
        else
            lmap
    })
    listmap <- IRanges::SplitDataFrameList(listmap)
    newmap <- IRanges::stack(listmap, "assay")
    newmap[["assay"]] <- factor(newmap[["assay"]])
    newmap
}
