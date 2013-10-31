corFinder <- function  # Calculate pair-wise correlations between samples using the expr() slots of a list of two ExpressionSets.
### This function acts as a wrapper around ComBat (sva package) and
### cor(), to calculate pairwise correlations within one or between
### two ExpressionSets.
(eset.pair,
### a list of ExpressionSets, with two elements.  If the two
### elements are identical, return the correlation matrix for pairs
### of samples in the first element.  If not identical, return pairs
### between the two elements.
separator=":",
### Separator between dataset name and sample name.  Dataset names are
### added to sample names to keep track of dataset of origin.
use.ComBat=TRUE,
### Use the sva::ComBat function for batch correction of the expr()
### data between the two datasets.
...
### Extra arguments passed to the cor() function.
 ){
    if((class(eset.pair) != "ExpressionSet")
       & (class(eset.pair) != "list" || length(eset.pair) > 2))
        stop("eset.pair should be a list of two esets")
    if(is.null(names(eset.pair)))
        names(eset.pair) <- paste("eset", 1:2, sep="")
    if( identical(class(eset.pair), "list") & !identical(eset.pair[[1]], eset.pair[[2]]) ){
        genes.intersect <- intersect(featureNames(eset.pair[[1]]), featureNames(eset.pair[[2]]))
        samples.intersect <- intersect(sampleNames(eset.pair[[1]]), sampleNames(eset.pair[[2]]))
        for (i in 1:length(eset.pair)){
            eset.pair[[i]] <- eset.pair[[i]][genes.intersect, samples.intersect]
            sampleNames(eset.pair[[i]]) <- paste(names(eset.pair)[i], sampleNames(eset.pair[[i]]), sep=separator)
        }
        ## Calculate correlation matrix for a pair of ExpressionSets:
        if(use.ComBat){
            big.matrix <- do.call(cbind, lapply(eset.pair, exprs))
            batch.var <- lapply(names(eset.pair), function(x) rep(x, ncol(eset.pair[[x]])))
            batch.var <- do.call(c, batch.var)
            big.matrix.combat <- sva::ComBat(big.matrix, mod=model.matrix(~(rep(1, length(batch.var)))), batch=batch.var)
            matrix.pair <- lapply(unique(batch.var), function(x) big.matrix.combat[, batch.var %in% x])
            names(matrix.pair) <- unique(batch.var)
        }else{
            matrix.pair <- lapply(eset.pair, exprs)
            names(matrix.pair) <- names(eset.pair)
        }
        cormat <- cor(matrix.pair[[1]], matrix.pair[[2]], ...)
    }else{
        ##Calculate correlation matrix for a single ExpressionSet:
        if(identical(class(eset.pair), "list")){
            matrix.one <- exprs(eset.pair[[1]])
            colnames(matrix.one) <- paste(names(eset.pair)[1], colnames(matrix.one), sep=separator)
        }else{
            matrix.one <- exprs(eset.pair)
        }
        cormat <- cor(matrix.one, ...)
        cormat[!upper.tri(cormat)] <- NA  ##NA for diagonal
    }
    return(cormat)
###   Returns a matrix of sample-wise Pearson Correlations.
}
