#' @param listmap A \code{list} class object containing names of either 
#' experiments, assays or features. 
#' @param type Any of the valid types of maps including colnames, rownames, 
#' and assays. 
#' @return A \code{DataFrame} class object of names
#' @describeIn mapToList Inverse of the listToMap function
#' @export listToMap
listToMap <- function(listmap, type = "colnames") {
  DFmap <- lapply(seq_along(listmap), FUN = function(i, x) {
    if (type == "colnames") {
      if (S4Vectors::isEmpty(x[i])) {
        S4Vectors::DataFrame(master = Rle(NA),
                             assay = NA,
                             assayname = Rle(names(x)[i]))
      } else {
        S4Vectors::DataFrame(master = Rle(x[[i]][, 1]),
                             assay = x[[i]][, 2],
                             assayname = Rle(names(x)[i]))
      }
    } else if (type == "rownames") {
      if (S4Vectors::isEmpty(x[i])) {
        S4Vectors::DataFrame(feature = NA,
                             assayname = Rle(names(x)[i]))
      } else {
        S4Vectors::DataFrame(feature = x[[i]][, 1],
                             assayname = Rle(names(x)[i]))
      }
    } else if (type == "assays") {
      S4Vectors::DataFrame(value = x[[i]],
                           assayname = names(x)[i])
    }
  }, x = listmap)
  newMap <- do.call(S4Vectors::rbind, DFmap)
  newMap <- newMap[!is.na(newMap[, 1]), ]
  return(newMap)
}
