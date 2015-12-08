#' plot features from two assays
#' 
#' @title plot features chosen from two assays in a MultiAssayExperiment
#' @description plot features chosen from two assays in a MultiAssayExperiment
#' @param fmla1 defines sources of y and x using names of Elist
#' @param fmla2 defines scatterplot using names of features in the Elist elements
#' @param mae the MultiAssayExperiment instance
#' @param ytx function transforming y before plot
#' @param xtx function transforming x before plot
#' @param \dots passed to plot
#' @export pwplot

pwplot = function(fmla1, fmla2, mae, ytx=force, xtx=force, ...) {
#
# use fmla1 with assays as components to identify 
#    two assays to regard as sources of y and x
# fmla2 indicates which features to plot
#
lf = as.list(fmla1)
nms = lapply(lf, as.character)
yel = Elist(mae)[[nms[[2]]]]
xel = Elist(mae)[[nms[[3]]]]
sy = samples(yel)
sx = samples(xel)
sb = intersect(sy,sx)
yel = yel[,sb]
xel = xel[,sb]
lf2 = lapply(as.list(fmla2), as.character)
ndf = data.frame( ytx(assay(yel)[ lf2[[2]], ]), xtx(assay(xel)[ lf2[[3]], ]) )
names(ndf) = c(lf2[[2]], lf2[[3]])
plot(fmla2, ndf, ...)
}
