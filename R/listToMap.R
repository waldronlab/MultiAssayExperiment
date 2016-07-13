#' @param listmap A \code{list} class object containing names of either 
#' experiments, assays or features. 
#' @param type Any of the valid types of maps including colnames, rownames, 
#' and assays. 
#' @return A \code{DataFrame} class object of names
#' @describeIn mapToList Inverse of the listToMap function
#' @export listToMap
listToMap <- function(listmap, type = "colnames") {
    if (!type %in% c("colnames", "rownames", "assays"))
        stop("Type not valid")
    DFmap <- lapply(seq_along(listmap), FUN = function(i, x) {
        if (type == "colnames") {
            if (S4Vectors::isEmpty(x[i])) {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     primary = NA,
                                     colname = NA)
            } else {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     primary = x[[i]][, 1],
                                     colname = x[[i]][, 2])
            }
        } else if (type == "rownames") {
            if (S4Vectors::isEmpty(x[i])) {
                S4Vectors::DataFrame(assay = factor(names(x)[i]),
                                     rowname = NA)
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
