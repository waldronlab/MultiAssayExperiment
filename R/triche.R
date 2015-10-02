# 
# #require(GenomicRanges) ## ensure that SummarizedExperiment has been declared
# setClassUnion("SummarizedExperimentOrNULL",c("SummarizedExperiment","NULL"))
# 
# ## slots: exprData, methData, geneData, exonData, lincData, mirData, cnvData
# setClass("TTMergedDataSet", contains=c("eSet"), # {{{
#          representation=representation(
#            exprData="SummarizedExperimentOrNULL",   # exprs, pvals?
#            methData="SummarizedExperimentOrNULL",   # betas, totals
#            geneData="SummarizedExperimentOrNULL",   # counts, RPKM, pvals?
#            exonData="SummarizedExperimentOrNULL",   # counts, RPKM
#            lincData="SummarizedExperimentOrNULL",   # counts, RPKM
#            mirData="SummarizedExperimentOrNULL",    # counts, RPKM
#            cnvData="SummarizedExperimentOrNULL"     # RLE CN, pvals?
#         )) # }}}
# 
# setValidity("TTMergedDataSet", function(object) { # {{{ check sample names
#   valid = TRUE
#   genome = NULL
#   samples = NULL
#   for( s in c('exprData', 'methData', 'geneData',
#               'exonData', 'lincData', 'mirData', 'cnvData') ) {
#     if( !empty(slot(object, s)) ) {
# 
#       # do all of the genomes match?
#       if(is(slot(object, s), 'SummarizedExperiment')) {
#         if(is.null(genome)) genome = unique(genome(slot(object, s)))
#         valid = valid && (genome == unique(genome(slot(object, s))))
#       }
# 
#       # do all of the samples match?
#       if(is(slot(object, s), 'ExpressionSet')) {
#         samples.s = sampleNames(slot(object, s))
#       } else if(is(slot(object, s), 'SummarizedExperiment')) {
#         samples.s = colnames(slot(object, s))
#       }
# 
#       if(is.null(samples)) samples = samples.s
#       valid = valid && (samples == samples.s)
# 
#     }
#   }
#   return( valid )
# }) # }}}
# 
# setMethod("show",signature(object="TTMergedDataSet"), function(object) { #
# callNextMethod()
#   available = c()
#   for( s in c('exprData', 'methData', 'geneData',
#               'exonData', 'lincData', 'mirData', 'cnvData') ) {
#     if(!empty(slot(object, s))) available = append(available, s)
#   }
#   if(empty(available)) cat('No experimental data has been added yet.\n')
#   else cat('Available data types:\n', paste(available, collapse=', '), '\n')
#   message("FIXME: add methods to automatically add columns of NAs to new data")
# })
