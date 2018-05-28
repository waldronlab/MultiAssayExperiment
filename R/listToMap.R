#' @param listmap A named \code{list} object containing \code{DataFrame}s
#' with "primary" and "colname" columns
#' @importFrom IRanges stack
#'
#' @return A \linkS4class{DataFrame} class object of names
#' @describeIn mapToList The inverse of the listToMap operation
#' @export listToMap
listToMap <- function(listmap) {
    if (is.null(names(listmap)))
        stop("'listmap' must be a named list")

    elementClass <- unique(vapply(listmap, class, character(1L)))
    if (!elementClass %in% c("DataFrame", "data.frame") ||
            length(elementClass) != 1L)
        stop("'listmap' elements are not all 'DataFrame' or 'data.frame'")

    if (elementClass == "data.frame")
        listmap <- lapply(listmap, S4Vectors::DataFrame)

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
