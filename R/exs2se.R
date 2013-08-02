exs2se = function(x, 
     assayname="exprs",
     fngetter = 
           function(z) rownames(exprs(z)),
#          function(z) fData(z)$probeset, 
     chipDbGetter = 
          function(z) { 
              clnanno = sub(".db", "", annotation(z))
              stopifnot(require(paste0(annotation(z), ".db"), character.only=TRUE) )
              get(paste0(annotation(z), ".db"))
              },
     probekeytype="PROBEID",
     duphandler=function(z) {
          if (any(isd <- duplicated(z[,"PROBEID"])))
              return(z[!isd,,drop=FALSE])
          z
          },
     signIsStrand = TRUE,
     ucsdChrnames = TRUE
 ) {
#
# first stab at ExpressionSet to SummarizedExperimentConverter
#
# acquire probe address source
 stopifnot(is(annopk <- chipDbGetter(x), "ChipDb"))
# acquire probe ids
 fn = fngetter(x)
# acquire chr, start and end
 locd = duphandler(fulls <- select(annopk, keytype=probekeytype, keys=fn,
           columns=c("CHR", "CHRLOC", "CHRLOCEND")))
 nfulls = na.omit(fulls)
 nmultiaddr = nrow(nfulls) - length(fn)
 rownames(locd) = locd[, probekeytype]
 locd = na.omit(locd)
 dropped = setdiff(fn, rownames(locd))
 if (length(dropped)>0) warning(paste("there were", length(dropped), 
    "addresses dropped owing to missing address information in bioc annotation"))
 locd = locd[intersect(rownames(locd), fn),]
 strand = rep("*", nrow(locd))
 if (signIsStrand) strand = ifelse(locd[,"CHRLOC"]>0, "+", "-")
 chpref = ""
 if (ucsdChrnames) chpref="chr"
 rowd = GRanges( paste0(chpref, locd[,"CHR"]), IRanges(
      abs(locd[,"CHRLOC"]), abs(locd[,"CHRLOCEND"])), strand=strand)
 names(rowd) = rownames(locd)
 metadata(rowd)$dropped = dropped
 metadata(rowd)$nmultiaddr = nmultiaddr
 hasrowd = match(names(rowd), fn, nomatch=0)
 ex = x[hasrowd,]
 stopifnot(nrow(exprs(ex)) == length(rowd))
# ad = SimpleList()
# ad[[assayname]] = exprs(ex)
 ed = SimpleList(initExptData=experimentData(ex))
 SummarizedExperiment(assays=SimpleList(exprs=exprs(ex)), 
       rowData=rowd, 
       colData=DataFrame(pData(ex)),
       exptData=ed
       )
}

setAs("ExpressionSet", "SummarizedExperiment",
   function(from) {
# coerce
     exs2se(from)
   })
