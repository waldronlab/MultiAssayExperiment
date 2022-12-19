#' @param listmap A named \code{list} object containing \code{DataFrame}s
#'   with "primary" and "colname" columns
#'
#' @param fill logical(1) Whether to fill the map with an empty \code{DataFrame}
#'   when empty elements are present in the input list
#'
#' @return A \linkS4class{DataFrame} class object of names
#' @describeIn mapToList The inverse of the listToMap operation
#' @export listToMap
listToMap <- function(listmap, fill = TRUE) {
    if (is.null(names(listmap)))
        stop("'listmap' must be a named list")
    if (!BiocBaseUtils::isTRUEorFALSE(fill))
        stop("'fill' must be a logical value")

    elementClass <- unique(vapply(listmap, class, character(1L)))
    elementNames <- names(listmap)

    if (length(elementClass) != 1L)
        stop("All 'listmap' elements must be of the same class")

    alldfs <- all(
        vapply(
            listmap,
            function(m) { is(m, "DataFrame") || is.data.frame(m) },
            logical(1L)
        )
    )
    if (!alldfs)
        stop("'listmap' elements are not all 'DataFrame' or 'data.frame'")

    if (!is(listmap, "SplitDataFrameList"))
        listmap <- IRanges::DataFrameList(listmap)

    listmap <- lapply(listmap, function(lmap) {
        if (isEmpty(lmap) && fill)
            DataFrame(primary = NA_character_, colname = NA_character_)
        else
            lmap
    })
    listmap <- as(listmap, "SplitDataFrameList")
    newmap <- IRanges::stack(listmap, "assay")
    newmap[["assay"]] <- factor(newmap[["assay"]], levels = elementNames)
    newmap
}
