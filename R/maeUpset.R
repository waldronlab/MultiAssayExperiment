#' create a generalized venn diagram analog for sample membership in multiple assays, using the upset algorithm in UpSetR
#'
#' @importFrom UpSetR upset
#' @param mae instance of \code{\link[MultiAssayExperiment]{MultiAssayExperiment-class}}
#' @param nsets integer number of sets to analyze, as in \code{\link[UpSetR]{upset}}
#' @param nintersects integer number of intersections to display
#' @param order character, either 'freq' or 'degree' to indicate order of intersection size display, as in \code{\link[UpSetR]{upset}}
#' @param idclip a function that operates on \code{colnames(mae)}, to remove potentially assay-specific token components; use \code{force} if no clipping is needed
#' @param \dots parameters passed to \code{\link[UpSetR]{upset}}
#' @examples
#' example(MultiAssayExperiment)
#' upset_samples(myMultiAssayExperiment, idclip=function(x) gsub("[a-z]", "", x))
#'
#' @export upset_samples
upset_samples = function(mae, nsets=length(experiments(mae)),
    nintersects=24, order="freq",
    idclip = function(x) substr(x, 1, 12), ... ) {
  maesn = colnames(mae)
  st = idclip(maesn[[1]])
  for (i in 2:length(maesn)) st = union(st, idclip(maesn[[i]]))
  nr = length(st)
  incid = matrix(0, nr=nr, nc=length(maesn))
  rownames(incid) = as.character(st)
  for (i in 1:length(maesn)) 
    incid[,i] = 1*(rownames(incid) %in% idclip(maesn[[i]]))
  colnames(incid) = names(experiments(mae))
  upset(data.frame(incid), nsets=nsets, nintersects=nintersects, 
      sets=colnames(incid), order=order, ...)
}
