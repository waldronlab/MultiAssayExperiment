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

    if (!elementClass %in% c("DataFrame", "data.frame") ||
            length(elementClass) != 1L)
        stop("'listmap' elements must all be 'DataFrame' or 'data.frame'")

    if (elementClass == "data.frame")
        listmap <- IRanges::SplitDataFrameList(
            lapply(listmap, S4Vectors::DataFrame))

    avector <- factor(rep(names(listmap), lengths(listmap)))
    cbind(assay = avector, unlist(listmap, use.names = FALSE))
}
