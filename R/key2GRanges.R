key2GRanges <- structure(function
### Convert an identifier (such as gene symbols) to a GRanges object, using
### a Bioconductor .db annotation database.
(x,
### character vector of the keys to convert
db="org.Hs.eg.db",  
### Name of the Bioconductor .db annotation package to use
keytype="SYMBOL",
### Must be a valid key type for the .db package, e.g. see keytypes(org.Hs.eg.db)
duphandler=function(z) {
  if (any(isd <- duplicated(z[, keytype])))
    return(z[!isd,,drop=FALSE])
  z},
### A function for handling multiple ranges returned for one key
  signIsStrand = TRUE,
### If TRUE, negative values are annotated as being on the - strand. 
  ucsdChrnames = TRUE
### If TRUE, use chr1 etc.
){
  ## acquire chr, start and end
  library(package=db, character.only=TRUE)
  locd = duphandler(fulls <- select(get(db), keytype=keytype, keys=x,
                                    columns=c("CHR", "CHRLOC", "CHRLOCEND")))
  nfulls = na.omit(fulls)
  nmultiaddr = nrow(nfulls) - length(x)
  rownames(locd) = locd[, keytype]
  locd = na.omit(locd)
  dropped = setdiff(x, rownames(locd))
  if (length(dropped)>0) warning(paste("there were", length(dropped), 
                                       "addresses dropped owing to missing address information in bioc annotation"))
  locd = locd[intersect(rownames(locd), x),]
  strand = rep("*", nrow(locd))
  if (signIsStrand) strand = ifelse(locd[,"CHRLOC"]>0, "+", "-")
  chpref = ""
  if (ucsdChrnames) chpref="chr"
  rowd = GRanges( paste0(chpref, locd[,"CHR"]), IRanges(
    abs(locd[,"CHRLOC"]), abs(locd[,"CHRLOCEND"])), strand=strand)
  names(rowd) = rownames(locd)
  metadata(rowd)$dropped = dropped
  metadata(rowd)$nmultiaddr = nmultiaddr
  rowd
### Returns a GRanges object
},ex=function(){
  key2GRanges(c("TP53", "BRCA1"))
  key2GRanges(c("1", "10"), keytype="ENTREZID")
})
