#' compute all pairwise linear models with optional transformations
#'
#' @param fmla A formula specifying assays using Elist element names
#' @param mae A MultiAssayExperiment instance
#' @param xtx a function that transforms the independent variable
#' @param ytx a function that transforms the dependent variable
#' @return a list with components 'mods' (lm results) and 'tslopes' (t stats for slope)
#' @export allLM_pw

allLM_pw = function(fmla, mae, xtx=force, ytx=force) {
#
# formula specifies dependent and independent assays
# form all regressions of ytx(dep) on xtx(indep) for all
# pairs of dependent and independent variables defined by
# assays for samples held in common
#
lf = as.list(fmla)
nms = lapply(lf, as.character)
yel = Elist(mae)[[nms[[2]]]]
xel = Elist(mae)[[nms[[3]]]]
sy = samples(yel)
sx = samples(xel)
sb = intersect(sy,sx)
yel = yel[,sb]
xel = xel[,sb]
vdf = as.matrix(expand.grid( names(rownames(yel)), 
    names(rownames(xel)), stringsAsFactors=FALSE ))
allf = apply(vdf, 1, function(x) as.formula(paste(x, collapse="~")))
alllm = foreach (i = 1:length(allf)) %dopar% {
  df = data.frame(ytx(assay(yel)[vdf[i,1],]), xtx(assay(xel)[vdf[i,2],]))
  names(df) = vdf[i,]
  lm(allf[[i]], data=df)
  }
names(alllm) = apply(vdf,1,function(x) paste(x, collapse="~"))
allts = lapply(alllm, function(x) summary(x)$coef[2, "t value"])
list(mods=alllm, tslopes=allts)
}

