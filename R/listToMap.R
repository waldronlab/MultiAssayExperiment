#' @param listmap A named \code{list} object containing either
#' experiments (\code{assay}s), samples (\code{colname}s) or
#' features (\code{rowname}s)
#' @param type Any of the valid types of maps including colnames, rownames,
#' and assays.
#' @return A \linkS4class{DataFrame} class object of names
#' @describeIn mapToList The inverse of the listToMap operation
#' @export listToMap
listToMap <- function(listmap, type = "colnames") {
    if (is.null(names(listmap)))
        stop("'listmap' must be a named list")
    type <- match.arg(type, c("colnames", "rownames", "assays"))
    DFmap <- lapply(seq_along(listmap), FUN = function(i, x) {
        if (type == "colnames") {
            if (S4Vectors::isEmpty(x[i])) {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     primary = NA_character_,
                                     colname = NA_character_)
            } else {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     primary = x[[i]][, 1],
                                     colname = x[[i]][, 2])
            }
        } else if (type == "rownames") {
            if (S4Vectors::isEmpty(x[i])) {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     rowname = NA_character_)
            } else {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     rowname = x[[i]])
            }
        } else if (type == "assays") {
            S4Vectors::DataFrame(assay = factor(names(x)[i]))
        }
    }, x = listmap)
    newMap <- do.call(S4Vectors::rbind, DFmap)
    newMap <- newMap[!is.na(newMap[, 1]), , drop = FALSE]
    return(newMap)
}
